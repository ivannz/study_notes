"""Function for testing"""
import numpy as np

from scipy.linalg import cholesky
from sklearn.metrics.pairwise import pairwise_kernels as kernel

from sklearn.utils import check_random_state

def rosenbrock(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `rosenbrock`"""
    val_ = 100 * np.sum((X[:, 1:] - X[:, :-1]**2)**2, axis=1)
    val_ += np.sum((X[:, :-1] - 1)**2, axis=1)
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def auckley(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `auckley`"""
    x, y = X[:, 0] * 4, X[:, 1] * 4
    val_ = -20 * np.exp(-0.2 * np.sqrt((x**2 + y**2) / 2))
    val_ -= np.exp((np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y)) / 2)
    val_ += 20 + np.e
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def eggholder(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `eggholder`"""
    x, y = X[:, 0] * 400, X[:, 1] * 400 + 47
    val_ = - y * np.sin(np.sqrt(np.abs(0.5*x + y))) - x * np.sin(np.sqrt(np.abs(x - y)))
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def levi(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `levi`"""
    x, y = X[:, 0] + 1, X[:, 1] + 1
    v1_ = 1 + np.sin(3 * np.pi * y)**2
    v2_ = 1 + np.sin(2 * np.pi * y)**2
    v1_ *= (x - 1)**2
    v2_ *= (y - 1)**2
    val_ = np.sin(3 * np.pi * x)**2 + v1_ + v2_
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def holder(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `holder`"""
    x, y = X[:, 0] * 4, X[:, 1] * 4
    val_ = - np.abs(np.sin(x) * np.cos(y) * np.exp(np.abs(1 - np.sqrt(x**2 + y**2) / np.pi)))
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def schaffer(X, random_state=None, scale=1.0, **kwargs):
    """Test function: `schaffer`"""
    x, y = 20 * (X[:, 0]**2), 20 * (X[:, 1]**2)
    val_ = 0.5 + (np.sin(x - y)**2 - 0.5) / (1 + (x + y) / 1000)
    if scale > 0:
        val_ /= val_.std()
        val_ *= scale
    return val_

def gaussian(X, size=None, scale=1.0, random_state=None, nugget=1e-6, metric="rbf", **kwargs):
    """Generate a realisation of a Gaussian process sampled at `X`.
    
    Parameters
    ----------
    nugget : float, optional, default=1e-6
        Cholesky decomposition is aaplicable to strictly positive definite
        matrices. Since Kernel matrices are positive-definite (semi-definite)
        `nugget` adds some regularization to the kernel, which is equvalent
        to mixing in some white noise with varaince `nugget` at each x in X.
    """
## Prepare the kernel matrix 
    Kxx = np.asfortranarray(kernel(X, metric=metric, **kwargs))
## Add some White noise and do the in-place cholesky decomposition.
    Kxx[np.diag_indices_from(Kxx)] += nugget
    cholesky(Kxx, lower=True, overwrite_a=True)
## Draw some independent gaussian variates.
    random_state = check_random_state(random_state)
    size_ = size if size is not None else tuple()
    return np.dot(Kxx, scale * random_state.normal(size=(Kxx.shape[1],) + size_))

def get_functions():
    """Returns a dictionary"""
    functions_ = [rosenbrock, auckley, eggholder,
                  levi, holder, schaffer, gaussian]
    return {fn.__name__: fn for fn in functions_}
