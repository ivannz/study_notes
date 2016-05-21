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
BASE_PATH = os.path.join(".", "exp_gauss_1d")
if not os.path.exists(BASE_PATH):
    os.mkdir(BASE_PATH)

n_jobs, verbose = -1, 0
parallel_ = Parallel(n_jobs=n_jobs, verbose=verbose)

## The levels and the random state
seeds_ = [0xB5066DBC, 0x98E8576F, 0x3161F88E, 0x08CCA9D9,]
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
grid_ = ParameterGrid(dict(size=[25, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600,],
                           nugget=true_nugget,
                           theta0=[1e+1, 1e+2, 1e+3, "auto"]))
## Initialize
kernel = 'rbf' # 'laplacian'
gp = GaussianProcess(beta0=0, normalize=False, corr='squared_exponential')

# Generate input
XX_test = np.linspace(0, 1, num=1001).reshape((-1, 1))
XX_train = random_state.uniform(size=(10000, 1))
XX = np.concatenate([XX_test, XX_train], axis=0)
test_ = np.s_[:XX_test.shape[0]]

experiment, batch_, dumps_ = list(), 1, list()
for noise_ in true_nugget:
    yy = gaussian(XX, scale=1.0, nugget=noise_, metric=kernel,
                  gamma=true_theta, random_state=random_state)
    if yy.ndim == 1:
        yy = yy.reshape((-1, 1))

    ## Split the pooled sample
    yy_train, yy_test = np.delete(yy, test_, axis=0), yy[test_].copy()
    for i_, par_ in enumerate(grid_):
        print i_, par_
        size_, nugget_, theta0_ = par_['size'], par_['nugget'], par_['theta0']

        tick_ = time.time()
        n_replications, replications = 25, list()
        while n_replications > 0:
            ## Draw random train sample
            train_ = random_state.choice(range(XX_train.shape[0]),
                                         size=size_, replace=False)
            X, y = XX_train[train_], yy_train[train_]

            ## Fit a GPR
            gp_ = clone(gp)
            gp_.nugget = nugget_
            if isinstance(theta0_, float):
                gp_.theta0 = theta0_
            elif theta0_ == "auto":
                gp_.thetaL, gp_.thetaU, gp_.theta0 = 1.0, 1e4, float(size_)
            gp_.fit(X, y)

            ## Compute the A, B matrices
            A, B, y_hat_, MM, loo_, A_loo, B_loo = \
                KRR_AB(X, y, XX_test, forecast=True,
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

            n_replications -= 1
            replications.append((yy_test[:, 0], y_hat_[:, 0],
                                 bp_width_, bp_hits_.mean(axis=0, keepdims=False),
                                 bf_width_, bf_hits_.mean(axis=0, keepdims=False),
                                 rrcm_width_, rrcm_hits_.mean(axis=0, keepdims=False),
                                 crr_width_, crr_hits_.mean(axis=0, keepdims=False),
                                 loo_rrcm_width_, loo_rrcm_hits_.mean(axis=0, keepdims=False),
                                 loo_crr_width_, loo_crr_hits_.mean(axis=0, keepdims=False)))
        tock_ = time.time()
        print "%0.3fsec"%(tock_-tick_,)

        key_ = "gaussian", noise_, theta0_, nugget_, size_
        result_ = tuple(np.stack([rep_[j] for rep_ in replications], axis=-1) for j in xrange(14))
        experiment.append((key_,) + result_)

        if len(experiment) >= 25:
            basename_ = os.path.join(BASE_PATH, "exp_gauss_1d%04d"%(batch_,))
            dumps_.append(_save(experiment, basename_, gz=9))
            experiment, batch_ = list(), batch_ + 1

if len(experiment) > 0:
    basename_ = os.path.join(BASE_PATH, "exp_gauss_1d%04d"%(batch_,))
    dumps_.append(_save(experiment, basename_, gz=9))
