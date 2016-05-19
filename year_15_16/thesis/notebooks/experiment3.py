"""None"""

import time
import os
import numpy as np

from sklearn.grid_search import ParameterGrid
from sklearn.base import clone

# from sklearn.utils import check_random_state
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcess

from scipy.stats import norm
from joblib import Parallel, delayed

from utils.state import _save
from utils.functions import gaussian

from utils.conformal import RRCM, CRR
from utils.KRR import KRR_AB

np.seterr(all="ignore")

def mkdirifnot(path):
    """Create a folder unless it exists.
    """
    if not os.path.exists(path):
        os.mkdir(path)
    return path

BASE_PATH = mkdirifnot(os.path.join(".", "exp3"))

## The levels and the random state
# random_state = np.random.RandomState(0x0BADA550)
# random_state = np.random.RandomState(0xB14FCA33)
# random_state = np.random.RandomState(0x1337BACE)
random_state = np.random.RandomState(0xCAFFE14E)

levels = np.asanyarray([0.01, 0.05, 0.10, 0.25])[::-1]


## Define the grid
grid_ = ParameterGrid(dict(dgp=[gaussian,],
                           size=[25, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600,],
                           # size=[150, 1500,],
                           nugget=[1e-6, 1e-2,],
                           theta0=[1e-1, 1, 1e+1, "auto"],
                           noise=[1e-6,]))

## Initialize
scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
gp = GaussianProcess(beta0=0, normalize=False, corr='squared_exponential')
kernel = 'rbf' # 'laplacian'

n_jobs, verbose = -1, 0
parallel_ = Parallel(n_jobs=n_jobs, verbose=verbose)

def _helper(y, A, B, proc=RRCM, levels=levels, parallel=None, n_jobs=1, verbose=0):
    if not isinstance(parallel, Parallel):
        parallel = Parallel(n_jobs=n_jobs, verbose=verbose)

## Construct the CKRR confidence interval: RRCM
    regions = parallel(delayed(proc)(A[k], B[k], levels=levels)
                       for k in xrange(y.shape[0]))

## See if the transformed test target valeus are with the conformal region
    hits_ = np.asarray(
        [[np.any(((int_[:, 0] <= target) & (target <= int_[:, 1]))).astype(float)
          for int_ in region]
         for target, region in zip(y, regions)])

    width_ = np.asarray(
        [[np.sum(int_[:, 1] - int_[:, 0]) for int_ in region] for region in regions])
    return hits_, width_


## Run: experiment #3
# nd = 2
# mesh_ = np.meshgrid(*nd*[np.linspace(-1, 1, num=51)])
# X_test = np.concatenate([ax_.reshape((-1, 1)) for ax_ in mesh_], axis=1)
X_test = np.linspace(0, 1, num=1001).reshape((-1, 1))
test_ = np.s_[:X_test.shape[0]]
## Draw f(x)
XX_train = random_state.uniform(size=(10000, 1))
XX = np.concatenate([X_test, XX_train], axis=0)
yy = gaussian(XX, scale=1.0, nugget=1e-6, metric=kernel,
              gamma=2.0, random_state=random_state)
if yy.ndim == 1:
    yy = yy.reshape((-1, 1))

## Split the pooled sample
yy_train, y_test = np.delete(yy, test_, axis=0), yy[test_].copy()
del XX, yy

experiment, batch_, dumps_ = list(), 1, list()
for i_, par_ in enumerate(grid_):
    print i_, par_
    n_replications, replications = 20, list()

    dgp_, size_, noise_ = par_['dgp'], par_['size'], par_['noise']
    nugget_, theta0_ = par_['nugget'], par_['theta0']

    tick_ = time.time()
    while n_replications > 0:
    ## START: one replication
    ## Draw random train sample
        train_ = random_state.choice(range(XX_train.shape[0]),
                                     size=size_, replace=False)
        X_train, y_train = XX_train[train_], yy_train[train_]

    ## Standardize the sample
        Xscl_, yscl_ = clone(scaler).fit(X_train), clone(scaler).fit(y_train)
        X_train_, X_test_ = Xscl_.transform(X_train), Xscl_.transform(X_test)
        y_train_, y_test_ = yscl_.transform(y_train), yscl_.transform(y_test)

    ## Fit a GPR
        gp_ = clone(gp)
        gp_.nugget = nugget_
        if isinstance(theta0_, float):
            gp_.theta0 = theta0_
        elif theta0_ == "auto":
            gp_.thetaL, gp_.thetaU, gp_.theta0 = 1e-4, 1e4, float(size_)
        gp_.fit(X_train_, y_train_)

    ## Compute the A, B matrices
        A, B, y_hat_, MM, loo_, A_loo, B_loo = KRR_AB(
            X_train_, y_train_, X_test_, forecast=True,
            nugget=gp_.nugget, metric=kernel, gamma=gp_.theta_[0])
        del loo_
    ## Inflate by the estimated magnitude
        MM *= gp_.sigma2

    ## Construct the Bayesian interval
        z_a = norm.ppf(1 - .5 * levels)
        half_width_ = np.sqrt(MM) * z_a[np.newaxis]
        b_bounds_ = yscl_.inverse_transform(
            np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1))
        b_width_ = b_bounds_[..., 1] - b_bounds_[..., 0]
        b_hits_ = ((b_bounds_[..., 0] <= y_test) & (y_test <= b_bounds_[..., 1])).astype(float)

    ## Construct the CKRR confidence interval: RRCM
        rrcm_hits_, rrcm_width_ = _helper(y_test_, A[0], B, proc=RRCM,
                                          levels=levels, parallel=parallel_)
        if yscl_.scale_ is not None:
            rrcm_width_ *= yscl_.scale_

    ## Construct the CKRR confidence interval: CCR-sided
        crr_hits_, crr_width_ = _helper(y_test_, A[0], B, proc=CRR,
                                        levels=levels, parallel=parallel_)
        if yscl_.scale_ is not None:
            crr_width_ *= yscl_.scale_

    ## Construct the CKRR confidence interval: RRCM
        loo_rrcm_hits_, loo_rrcm_width_ = _helper(y_test_, A_loo[0], B_loo, proc=RRCM,
                                                  levels=levels, parallel=parallel_)
        if yscl_.scale_ is not None:
            loo_rrcm_width_ *= yscl_.scale_

    ## Construct the CKRR confidence interval: CCR-sided
        loo_crr_hits_, loo_crr_width_ = _helper(y_test_, A_loo[0], B_loo, proc=CRR,
                                                levels=levels, parallel=parallel_)
        if yscl_.scale_ is not None:
            loo_crr_width_ *= yscl_.scale_

    #     rrcm_regions = parallel_(delayed(RRCM)(A[0, k], B[k], levels=levels)
    #                              for k in xrange(y_test.shape[0]))

    # ## See if the transformed test target valeus are with the conformal region
    #     rrcm_hits_ = np.asarray(
    #         [[np.any(((int_[:, 0] <= y) & (y <= int_[:, 1]))).astype(float)
    #           for int_ in region]
    #          for y, region in zip(y_test_, rrcm_regions)])

    #     rrcm_width_ = np.asarray(
    #         [[np.sum(int_[:, 1] - int_[:, 0]) for int_ in region] for region in rrcm_regions])
    #     if yscl_.scale_ is not None:
    #         rrcm_width_ *= yscl_.scale_

    ## END: one replication
        n_replications -= 1
    
        y_test_hat_ = yscl_.inverse_transform(y_hat_)
        replications.append((y_test[:, 0], y_test_hat_[:, 0],
                             b_width_, b_hits_.mean(axis=0, keepdims=False),
                             rrcm_width_, rrcm_hits_.mean(axis=0, keepdims=False),
                             crr_width_, crr_hits_.mean(axis=0, keepdims=False),
                             loo_rrcm_width_, loo_rrcm_hits_.mean(axis=0, keepdims=False),
                             loo_crr_width_, loo_crr_hits_.mean(axis=0, keepdims=False)))
    tock_ = time.time()
    print "%0.3fsec"%(tock_-tick_,)

## Consolidate the simultions
    # b_cov_, rrcm_cov_, crr_cov_, loo_rrcm_cov_, loo_crr_cov_ = \
    #     [np.concatenate([rep_[j] for rep_ in replications], axis=0) for j in xrange(5)]
    # key_ = dgp_.__name__, noise_, theta0_, nugget_, size_
    # experiment.append((key_, b_cov_, rrcm_cov_, crr_cov_, loo_rrcm_cov_, loo_crr_cov_))

    key_ = dgp_.__name__, noise_, theta0_, nugget_, size_
    result_ = tuple(np.stack([rep_[j] for rep_ in replications], axis=-1) for j in xrange(12))
    experiment.append((key_,) + result_)

    if len(experiment) >= 25:
        basename_ = os.path.join(BASE_PATH, "exp3_%04d"%(batch_,))
        dumps_.append(_save(experiment, basename_, gz=9))
        experiment, batch_ = list(), batch_ + 1

if len(experiment) > 0:
    basename_ = os.path.join(BASE_PATH, "exp3_%04d"%(batch_,))
    dumps_.append(_save(experiment, basename_, gz=9))

# y_test, y_test_hat_, b_cov_, b_width_, rrcm_cov_, rrcm_width_, crr_cov_, \
# crr_width_, loo_rrcm_cov_, loo_rrcm_width_, loo_crr_cov_, loo_crr_width_ = \