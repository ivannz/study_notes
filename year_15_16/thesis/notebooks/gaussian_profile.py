"""None"""

import time
import os

import numpy as np

from sklearn.grid_search import ParameterGrid
from sklearn.base import clone

from sklearn.gaussian_process import GaussianProcess

from scipy.stats import norm
from joblib import Parallel, delayed

from utils.state import _save
from utils.functions import gaussian

from utils.conformal import RRCM, CRR
from utils.KRR import KRR_AB

np.seterr(all="ignore")
BASE_PATH = os.path.join(".", "prof_gauss")
if not os.path.exists(BASE_PATH):
    os.mkdir(BASE_PATH)

n_jobs, verbose = -1, 0
parallel_ = Parallel(n_jobs=n_jobs, verbose=verbose)

## The levels and the random state
seeds_ = [0x4A04E61B, 0x7A5F2F22, 0x52B4A070, 0x470BB766,]
random_state = np.random.RandomState(seeds_[1])

levels = np.asanyarray([0.01, 0.05, 0.10, 0.25])[::-1]

## helpers
def _helper(y, A, B, proc=RRCM, levels=levels, parallel=None, n_jobs=1, verbose=0):
    if not isinstance(parallel, Parallel):
        parallel = Parallel(n_jobs=n_jobs, verbose=verbose)

    regions = parallel(delayed(proc)(A[k], B[k], levels=levels)
                       for k in xrange(y.shape[0]))

    hits_ = np.asarray(
        [[np.any(((int_[:, 0] <= target) & (target <= int_[:, 1]))).astype(float)
          for int_ in region]
         for target, region in zip(y, regions)])

    width_ = np.asarray(
        [[np.sum(int_[:, 1] - int_[:, 0]) for int_ in region] for region in regions])
    
    bounds_ = np.asarray(
        [[[int_[:, 0].min(), int_[:, 1].max()] for int_ in region] for region in regions])

    return hits_, width_, bounds_

## Define the grid
true_theta = 100.0
true_nugget = [1e-6, 1e-1,]
sizes = [50, 200,]

grid_ = ParameterGrid(dict(nugget=true_nugget,
                           theta0=[1e+2, "auto"]))
## Initialize
kernel = 'rbf'
gp = GaussianProcess(beta0=0, normalize=False, corr='squared_exponential')
# kernel = 'laplacian'
# gp = GaussianProcess(beta0=0, normalize=False, corr='absolute_exponential')

# Generate input
XX_test = np.linspace(0, 1, num=501).reshape((-1, 1))
test_ = np.s_[:XX_test.shape[0]]

experiment, batch_, dumps_ = list(), 1, list()
for size_ in sizes:
    XX_train = np.linspace(0.05, 0.95, num=size_ + 1).reshape((-1, 1))
    XX = np.concatenate([XX_test, XX_train], axis=0)

    for noise_ in true_nugget:
        yy = gaussian(XX, scale=1.0, nugget=noise_, metric=kernel,
                      gamma=true_theta, random_state=random_state)
        if yy.ndim == 1:
            yy = yy.reshape((-1, 1))

        ## Split the pooled sample
        yy_train, yy_test = np.delete(yy, test_, axis=0), yy[test_].copy()

        for i_, par_ in enumerate(grid_):
            nugget_, theta0_ = par_['nugget'], par_['theta0']
            key_ = "gaussian", noise_, theta0_, nugget_, size_
            print i_, key_

            tick_ = time.time()
            ## Fit a GPR
            gp_ = clone(gp)
            gp_.nugget = nugget_
            if isinstance(theta0_, float):
                gp_.theta0 = theta0_
            elif theta0_ == "auto":
                gp_.thetaL, gp_.thetaU, gp_.theta0 = 1.0, 1e4, float(size_)
            gp_.fit(XX_train, yy_train)

            ## Compute the A, B matrices
            A, B, y_hat_, MM, loo_, A_loo, B_loo = \
                KRR_AB(XX_train, yy_train, XX_test, forecast=True,
                       nugget=gp_.nugget, metric=kernel, gamma=gp_.theta_[0])
            del loo_

            ## Construct the CKRR confidence interval: RRCM
            rrcm_hits_, rrcm_width_, rrcm_bounds_ = \
                _helper(yy_test, A[0], B, proc=RRCM,
                        levels=levels, parallel=parallel_)

            ## Construct the CKRR confidence interval: CCR-sided
            crr_hits_, crr_width_, crr_bounds_ = \
                _helper(yy_test, A[0], B, proc=CRR,
                        levels=levels, parallel=parallel_)

            ## Construct the CKRR confidence interval: RRCM
            loo_rrcm_hits_, loo_rrcm_width_, loo_rrcm_bounds_ = \
                _helper(yy_test, A_loo[0], B_loo, proc=RRCM,
                        levels=levels, parallel=parallel_)

            ## Construct the CKRR confidence interval: CCR-sided
            loo_crr_hits_, loo_crr_width_, loo_crr_bounds_ = \
                _helper(yy_test, A_loo[0], B_loo, proc=CRR,
                        levels=levels, parallel=parallel_)

            ## Construct the GPR forecast interval
            z_a = norm.ppf(1 - .5 * levels)
            half_width_ = np.sqrt(MM * gp_.sigma2) * z_a[np.newaxis]
            bf_bounds_ = np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1)
            bf_width_ = bf_bounds_[..., 1] - bf_bounds_[..., 0]
            bf_hits_ = ((bf_bounds_[..., 0] <= yy_test)
                        & (yy_test <= bf_bounds_[..., 1])).astype(float)

            ## Construct the GPR prediction interval
            half_width_ = np.sqrt((MM - gp_.nugget) * gp_.sigma2) * z_a[np.newaxis]
            bp_bounds_ = np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1)
            bp_width_ = bp_bounds_[..., 1] - bp_bounds_[..., 0]
            bp_hits_ = ((bp_bounds_[..., 0] <= yy_test)
                        & (yy_test <= bp_bounds_[..., 1])).astype(float)

            replication = (yy_test[:, 0], y_hat_[:, 0],
                           bp_bounds_, bf_bounds_,
                           rrcm_bounds_, crr_bounds_,
                           loo_rrcm_bounds_, loo_crr_bounds_,)

            tock_ = time.time()
            print "%0.3fsec"%(tock_-tick_,)

            result_ = tuple(rep_ for rep_ in replication)
            experiment.append((key_,) + result_)

basename_ = os.path.join(BASE_PATH, "prof_gauss%04d"%(batch_,))
dumps_.append(_save(experiment, basename_, gz=9))
