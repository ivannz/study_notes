"""None"""
import numpy as np

from scipy.stats import norm
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.utils import check_random_state
from sklearn.gaussian_process import GaussianProcess
from sklearn.preprocessing import StandardScaler

from .KRR import KRR_AB

def run_ckrr_mc_experiment(gdp, levels, ccr_proc, nd=1, ng=1001,
                           n_replications=1, size=100,
                           nugget=1e-6, theta0=1e-1, use_loo=False,
                           random_state=None):
    """One experiment. [0, 1]^d
    """
## Use RBF
    kernel = 'rbf'
    gp = GaussianProcess(beta0=0, theta0=theta0,
                         normalize=False, nugget=nugget,
                         corr='squared_exponential')

    scl = StandardScaler(copy=True, with_mean=True, with_std=True)

## Fix a test sample (for comparability)
    mesh_ = np.meshgrid(*nd*[np.linspace(0, 1, num=ng)])
    X_test = np.concatenate([ax_.reshape((-1, 1)) for ax_ in mesh_], axis=1)
    test_ = np.s_[:X_test.shape[0]]

## The `levels` quantiles of the standard normal
    z_a = norm.ppf(1 - .5 * levels)

## Initialize the parallel backend
    parallel_ = Parallel(n_jobs=-1, verbose=1)
    c_proc = delayed(ccr_proc)

    replications = list()
    random_state = check_random_state(random_state)
    while n_replications > 0:
## START: one replication
## Generate a train input sample
        X_train = random_state.uniform(size=(size, nd))

## Pool the inputs and produce targets
        XX = np.concatenate([X_test, X_train], axis=0)
        yy = gdp(XX)
        if yy.ndim == 1:
            yy = yy.reshape((-1, 1))

## Split the pooled sample
        y_train, y_test = np.delete(yy, test_, axis=0), yy[test_]

## Learning setting: standardize the datasets
        X_scl_ = clone(scl).fit(X_train)
        y_scl_ = clone(scl).fit(y_train)

        X_train_, X_test_ = X_scl_.transform(X_train), X_scl_.transform(X_test)
        y_train_, y_test_ = y_scl_.transform(y_train), y_scl_.transform(y_test)

## Fit a GPR
        gp_ = clone(gp).fit(X_train_, y_train_)

## Compute the A-B vectors
        A, B, y_hat_, MM, loo_ = KRR_AB(X_train_, y_train_, X_test_, loo=use_loo,
                                        forecast=True, nugget=gp_.nugget,
                                        metric=kernel, gamma=gp_.theta_[0])
        del loo_
## Inflate by the estimated magnitude
        MM *= gp_.sigma2

## Construct the Bayesian interval
        half_width_ = np.sqrt(MM) * z_a[np.newaxis]
        b_bounds_ = y_scl_.inverse_transform(
            np.stack([y_hat_ - half_width_, y_hat_ + half_width_], axis=-1))
        b_width_ = b_bounds_[..., 1] - b_bounds_[..., 0]
        b_hits_ = ((b_bounds_[..., 0] <= y_test) & (y_test <= b_bounds_[..., 1])).astype(float)

## Construct the CKRR confidence interval
        regions = parallel_(c_proc(A[0, k], B[k], levels=levels)
                            for k in xrange(y_test.shape[0]))

## return the convex hull of the confidence region (this increases its coverage a bit)
        # c_bounds_ = y_scl_.inverse_transform(np.asarray(
        #     [[[int_[0, 0], int_[-1, -1]] for int_ in region] for region in regions]))

## See if the transformed test target valeus are with the conformal region
        c_hits_ = np.asarray(
            [[np.any(((int_[:, 0] <= y) & (y <= int_[:, 1]))).astype(float)
              for int_ in region]
             for y, region in zip(y_test_, regions)])

        c_width_ = np.asarray(
            [[np.sum(int_[:, 1] - int_[:, 0]) for int_ in region] for region in regions])
        if y_scl_.scale_ is not None:
            c_width_ *= y_scl_.scale_

## END: one replication
        n_replications -= 1

## Compute the abs_error
        y_test_hat_ = y_scl_.inverse_transform(y_hat_)
        replications.append((y_test, y_test_hat_, b_width_, c_width_,
                             b_hits_.mean(axis=0, keepdims=True),
                             c_hits_.mean(axis=0, keepdims=True)))
    return X_test, replications
