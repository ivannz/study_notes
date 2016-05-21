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
from utils.functions_1d import f6, heaviside, pressure2

from utils.conformal import RRCM, CRR
from utils.KRR import KRR_AB

np.seterr(all="ignore")
BASE_PATH = os.path.join(".", "prof_nongauss")
if not os.path.exists(BASE_PATH):
    os.mkdir(BASE_PATH)

n_jobs, verbose = -1, 0
parallel_ = Parallel(n_jobs=n_jobs, verbose=verbose)

## The levels and the random state
seeds_ = [0x4A04E61B, 0x7A5F2F22, 0x52B4A070, 0x470BB766,]
random_state = np.random.RandomState(seeds_[0])

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

true_nugget = [1e-6, 1e-1,]
test_functions = [f6, heaviside, pressure2]

grid_ = ParameterGrid(dict(func=test_functions,
                           size=[200, 1200,],
                           nugget=true_nugget,
                           noise=true_nugget,
                           theta0=[1e+2, "auto"]))
## Initialize
kernel = 'rbf' # 'laplacian'
gp = GaussianProcess(beta0=0, normalize=False, corr='squared_exponential')

# Generate input
XX_test = np.linspace(0, 1, num=501).reshape((-1, 1))
test_ = np.s_[:XX_test.shape[0]]

experiment, batch_, dumps_ = list(), 1, list()
for i_, par_ in enumerate(grid_):
    print i_, par_
    size_, nugget_, theta0_ = par_['size'], par_['nugget'], par_['theta0']
    func_, noise_ = par_['func'], par_['noise']

    tick_ = time.time()
    n_replications, replications = 1, list()
    while n_replications > 0:
        ## Draw random train sample
        X = random_state.uniform(size=(size_, 1))

        ## Draw f(x)
        XX = np.concatenate([XX_test, X], axis=0)
        yy = func_(XX)
        if yy.ndim == 1:
            yy = yy.reshape((-1, 1))

        if noise_ > 0:
            yy += random_state.normal(size=yy.shape) * noise_

        ## Split the pooled sample
        y, y_test = np.delete(yy, test_, axis=0), yy[test_]

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
            _helper(y_test, A[0], B, proc=RRCM,
                    levels=levels, parallel=parallel_)

        ## Construct the CKRR confidence interval: CCR-sided
        crr_hits_, crr_width_, crr_bounds_ = \
            _helper(y_test, A[0], B, proc=CRR,
                    levels=levels, parallel=parallel_)

        ## Construct the CKRR confidence interval: RRCM
        loo_rrcm_hits_, loo_rrcm_width_, loo_rrcm_bounds_ = \
            _helper(y_test, A_loo[0], B_loo, proc=RRCM,
                    levels=levels, parallel=parallel_)

        ## Construct the CKRR confidence interval: CCR-sided
        loo_crr_hits_, loo_crr_width_, loo_crr_bounds_ = \
            _helper(y_test, A_loo[0], B_loo, proc=CRR,
                    levels=levels, parallel=parallel_)

        ## Construct the GPR forecast interval
        z_a = norm.ppf(1 - .5 * levels)
        half_width_ = np.sqrt(MM * gp_.sigma2) * z_a[np.newaxis]
        bf_bounds_ = np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1)
        bf_width_ = bf_bounds_[..., 1] - bf_bounds_[..., 0]
        bf_hits_ = ((bf_bounds_[..., 0] <= y_test)
                    & (y_test <= bf_bounds_[..., 1])).astype(float)

        ## Construct the GPR prediction interval
        half_width_ = np.sqrt((MM - gp_.nugget) * gp_.sigma2) * z_a[np.newaxis]
        bp_bounds_ = np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1)
        bp_width_ = bp_bounds_[..., 1] - bp_bounds_[..., 0]
        bp_hits_ = ((bp_bounds_[..., 0] <= y_test)
                    & (y_test <= bp_bounds_[..., 1])).astype(float)

        n_replications -= 1
        replications.append((y_test[:, 0], y_hat_[:, 0],
                             bp_bounds_, bf_bounds_,
                             rrcm_bounds_, crr_bounds_,
                             loo_rrcm_bounds_, loo_crr_bounds_,))

    tock_ = time.time()
    print "%0.3fsec"%(tock_-tick_,)

    key_ = func_.__name__, noise_, theta0_, nugget_, size_
    result_ = tuple(np.stack([rep_[j] for rep_ in replications], axis=-1) for j in xrange(8))
    experiment.append((key_,) + result_)

basename_ = os.path.join(BASE_PATH, "prof_nongauss_%04d"%(batch_,))
dumps_.append(_save(experiment, basename_, gz=9))
