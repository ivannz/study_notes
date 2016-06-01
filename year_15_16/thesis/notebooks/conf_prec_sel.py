"""None"""
import os

import numpy as np
import pandas as pd

from scipy.linalg import inv
from sklearn.metrics import pairwise_kernels
from sklearn.grid_search import ParameterGrid

from scipy.stats import norm
from joblib import Parallel, delayed

from utils.functions_1d import triangle

from utils.conformal import RRCM

from utils.state import _save

np.seterr(all="ignore")
BASE_PATH = os.path.join(".", "conf_prec_sel")
if not os.path.exists(BASE_PATH):
    os.mkdir(BASE_PATH)

levels = np.asanyarray([0.01, 0.05, 0.10, 0.25])[::-1]
levels_ = pd.Index(["%2.1f"%(100*lvl,) for lvl in levels], name="alpha")

random_state = np.random.RandomState(0x4A04B766)
metric = 'rbf'

master_grid = ParameterGrid(dict(noise=[1e-6, 5e-2, 5e-1,],
                                 size=[101, 251, 501,]))

nugget_grid = [1e-6, 1e-1,]
theta_grid = np.logspace(-2.0, 5.0, num=71)
z_a = norm.ppf(1 - .5 * levels)

n_replications = 1
experiment, batch_, dumps_ = list(), 1, list()
for par_ in master_grid:
    noise_, size_ = par_["noise"], par_["size"]

    X = np.linspace(-.5, 1.5, num=size_).reshape((-1, 1))
    y = triangle(X)
    y = y + random_state.normal(size=(y.shape[0], n_replications)) * noise_

# Compute the GPR and conformal loo-interval widths
    parallel_ = Parallel(n_jobs=-1, verbose=10)
    for loo_ in [True, False,]:
        for nugget_ in nugget_grid:
            loo_sigma_, y_hat_, A, B = list(), list(), list(), list()
            for theta_ in theta_grid:
                Kxx_ = np.asfortranarray(pairwise_kernels(X, metric=metric, gamma=theta_))
                Kxx_.flat[::Kxx_.shape[0] + 1] += nugget_

                inv(Kxx_, overwrite_a=True, check_finite=False)

                Mi_ = np.diag(Kxx_).reshape((1, -1))
                sigma2_ = np.diag(np.dot(y.T, np.dot(Kxx_, y))) / Kxx_.shape[0]
                loo_sigma2_ = sigma2_.reshape((-1, 1)) / Mi_
                if loo_: # leave-one-out predictions and residuals
                    yQM_ = np.dot(y.T, Kxx_) / Mi_
                    B_ = Kxx_ / Mi_
                    hat_ = y.T - yQM_
                else:    # in-sample predictions and residuals
                    yQM_ = np.dot(y.T, Kxx_)
                    hat_ = y.T - yQM_ / Mi_
                    yQM_ *= nugget_
                    B_ = nugget_ * Kxx_
                # Enqueue
                A.append(yQM_[:, np.newaxis] - B_[np.newaxis] * y[np.newaxis].T)
                loo_sigma_.append(np.sqrt(loo_sigma2_))
                y_hat_.append(hat_)
                B.append(B_)

            # Collect
            A, B = np.stack(A, axis=1), np.stack(B, axis=0)
            loo_sigma_ = np.stack(loo_sigma_, axis=1)
            y_hat_ = np.stack(y_hat_, axis=1)

            ## Run the conformal procedure on each column of A-B
            results_ = parallel_(delayed(RRCM)(A[r, t, n], B[t, n], levels=levels, n=n)
                                 for r in xrange(y.shape[1])
                                 for t in xrange(len(theta_grid))
                                 for n in xrange(X.shape[0]))

            rrcm_intervals_ = np.reshape(np.array([(res_[0, 0], res_[-1, -1])
                                                   for result_ in results_ for res_ in result_]),
                                         (y.shape[1], len(theta_grid), X.shape[0], len(levels), 2))

            gpr_hwidth_ = loo_sigma_[..., np.newaxis] * z_a.reshape((1, 1, 1, -1))
            gpr_intervals_ = np.stack([y_hat_[..., np.newaxis] - gpr_hwidth_,
                                       y_hat_[..., np.newaxis] + gpr_hwidth_], axis=-1)

            experiment.append((noise_, size_, nugget_, loo_,
                               (X, y, y_hat_, gpr_intervals_, rrcm_intervals_)))


            if len(experiment) >= 25:
                basename_ = os.path.join(BASE_PATH, "exp_triangle_1d%04d"%(batch_,))
                dumps_.append(_save(experiment, basename_, gz=9))
                experiment, batch_ = list(), batch_ + 1

if len(experiment) > 0:
    basename_ = os.path.join(BASE_PATH, "exp_triangle_1d%04d"%(batch_,))
    dumps_.append(_save(experiment, basename_, gz=9))

print dumps_

