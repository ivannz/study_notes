"""Procedures for Kernel Ridge Regression.

Since only $A$ and $B$ vector are required for construction of Conformal Confidence
Region for Regression construction, we prepare a special routine to compute them in bulk:
$$ B'e_i
    = \\begin{pmatrix}
        - Q_X K_{Xz_i} \\\\
        1
    \\end{pmatrix} m_i^{-1} a \\,, $$
and
$$ A'e_i
    = \\begin{pmatrix}
        Q_X y + Q_X K_{Xz_i} m_i^{-1} K_{z_iX} Q_X y \\\\
        - m_i^{-1} K_{z_iX} Q_X y
    \\end{pmatrix} a
    = \\begin{pmatrix} Q_X y \\\\ 0 \\end{pmatrix} a
    - \\begin{pmatrix}
        - Q_X K_{Xz_i} \\\\
        1
    \\end{pmatrix} m_i^{-1} K_{z_iX} Q_X y a
    = \\begin{pmatrix} a Q_X y \\\\ 0 \\end{pmatrix}
    - B'e_i K_{z_iX} Q_X y \\,, $$
with
$$ m_i = a + K(z_i, z_i) - K_{z_iX} Q_X K_{Xz_i} \\,. $$

Note that $a Q_x = I_n - K_x Q_x$.

LOO residuals are computed using the following result: for all $i=1,\\ldots, n$ it is
true that $\\hat{r}_i = e_i' a Q_X y$, which is given by 
$$ \\hat{r}_i
    = a m_i^{-1} \\bigl(y_i - k_{-i}(x_i)Q_{-i}y_{-i} \\bigr)
     = a m_i^{-1} \\hat{r}_{i\\vert -i}
    \\,, $$
using the block inversion of a row-columns permuted matrix $Q_X$. In a compacter matrix
form this is given by
$$ \\hat{r} = a \\mathop{\\text{diag}}(Q_X) \\hat{r}_{\\text{loo}} \\,, $$
which, when all $m_i$ are non-zero, is equivalent to;
$$ \\hat{r}_{\\text{loo}}
    = a^{-1} \\mathop{\\text{diag}}(Q_X)^{-1} \\hat{r}
    = a^{-1} \\mathop{\\text{diag}}(Q_X)^{-1} a Q_X y
    = \\mathop{\\text{diag}}(Q_X)^{-1} Q_X y
    \\,. $$
"""
import numpy as np

from scipy.linalg import cholesky, solve_triangular
from scipy.linalg.lapack import dtrtri

from sklearn.metrics.pairwise import pairwise_kernels

def KRR_AB(X, y, Z, forecast=True, nugget=1.0, metric="rbf", **kwargs):
    """Compute the A, B vectors sufficient for regression CCR construction.
    A, B have shape [m, n+1].
    X [N, ...], y [N, M], Z [K, ...]
    A [M, K, N+1]
    B [K, N+1]
    Returns
    -------
    A, ndarray [M, K, N+1]
        A[i, j, :] is the vector 
    """
    if Z.shape[0] < 1:
        raise ValueError("""At least one test object is required.""")
## Normalize the target vector
    if y.ndim < 2:
        y = y.reshape((-1, 1))
## Obtain AB vectors for the regression CCR problem
    Kxx = np.asfortranarray(pairwise_kernels(X, metric=metric, **kwargs))
    Kxx[np.diag_indices_from(Kxx)] += nugget
## Get the in-place lower triangular Cholesky decomposition (in-place works if
##  the array is in fortran layout).
    ## K_{xx} + a I_x = C C' ; K_{zx} Q_X y = (C^{-1} K_{xz})' (C^{-1} y)
    cholesky(Kxx, lower=True, overwrite_a=True)
    Ci_y = solve_triangular(Kxx, y.copy('F'), lower=True, overwrite_b=True)
## Compute the border vector and the corner value
    Kxz = np.asfortranarray(pairwise_kernels(X, Z, metric=metric, **kwargs))
    Kxz = solve_triangular(Kxx, Kxz, lower=True, overwrite_b=True)
## Compute the KzxQx product and diagonalize M
    if forecast:
        Kzz = np.diag(pairwise_kernels(Z, metric=metric, **kwargs)) + nugget
    else:
        Kzz = np.diag(pairwise_kernels(Z, metric=metric, **kwargs))
    neg_m = - Kzz + np.einsum("ji,ji->i", Kxz, Kxz)
    neg_m = neg_m[np.newaxis]
    del Kzz
## Start preparing the A matrix:
    y_z = np.dot(Ci_y.T, Kxz)
    solve_triangular(Kxx, Kxz, lower=True, trans=1, overwrite_b=True)
    solve_triangular(Kxx, Ci_y, lower=True, trans=1, overwrite_b=True)
    # Ci_y = Q_X y; Kxz = Q_X K_XZ
## Invert the Cholseky matrix and get the diagonal of L^{-1}' L^{-1}
    dtrtri(Kxx, overwrite_c=1, lower=1)
    diag_ = np.einsum("ji,ji->i", Kxx, Kxx, dtype=np.float128)
    del Kxx
## Now extract the B vector
    B = np.concatenate([Kxz / neg_m, (-1.0) / neg_m], axis=0)
## Now it's A'turn: column-wise multiplication of B
    A = (- B[np.newaxis]) * y_z[:, np.newaxis]
    A[..., :-1, :] += Ci_y.T[..., np.newaxis]
    Ci_y /= diag_[:, np.newaxis]

    A, B = A.swapaxes(-2, -1), B.T
    A_loo, B_loo = A.copy(), B.copy()

    A *= nugget
    B *= nugget

    # Kxz = Q_X K_XZ
    Kxz *= -B_loo[:, :-1].T ## NxK
    # Kxz = (Q_X K_XZ)^2 / M_Z
    Kxz += diag_[:, np.newaxis]
## `diag` contains the reciprocals of the leverages (up to `nugget` scaling)
    A_loo[..., :-1] /= Kxz.T
    B_loo[..., :-1] /= Kxz.T
    A_loo[..., -1] *= -neg_m[0]
    B_loo[..., -1] *= -neg_m[0]
    del Kxz
## return A, B, KzxQx, M and LOO-residuals: only M depends on sigma2.
    return A, B, y_z.T, - neg_m.T, Ci_y, A_loo, B_loo

def KRR_loo(X, y, nugget=1.0, metric="rbf", **kwargs):
    """A dedicated procedure to compute the LOO-residuals over
    the train dataset (X, y).
    """
## Normalize the target vector
    if y.ndim < 2:
        y = y.reshape((-1, 1))
## Obtain AB vectors for the regression CCR problem
    Kxx = np.asfortranarray(pairwise_kernels(X, metric=metric, **kwargs))
    Kxx[np.diag_indices_from(Kxx)] += nugget
## Get the in-place lower triangular Cholesky decomposition (in-place works if
##  the array is in fortran layout).
    ## K_{xx} + aI_x = C C' ; K_{zx} Q_X y = (C^{-1} K_{xz})' (C^{-1} y)
    cholesky(Kxx, lower=True, overwrite_a=True)
    Ci_y = solve_triangular(Kxx, y.copy('F'), lower=True, overwrite_b=True)
## Invert the Cholseky matrix and get the diagonal of L^{-1}' L^{-1}
    dtrtri(Kxx, overwrite_c=1, lower=1)
    diag_ = np.einsum("ji,ji->i", Kxx, Kxx, dtype=np.float128)
    return np.dot(Kxx.T, Ci_y) / diag_[:, np.newaxis]


def get_AB(X, y, Z, kernel, nugget, use_loo=False, **kwargs):
## Normalize the target vector
    if y.ndim < 2:
        y = y.reshape((-1, 1))
## Obtain AB vectors for the regression CCR problem
    Kxx = np.asfortranarray(kernel(X, **kwargs))
    Kxx.flat[::Kxx.shape[0]+1] += nugget
## Get the in-place lower triangular Cholesky decomposition (in-place works if
##  the array is in fortran layout).
    ## K_{xx} + a I_x = C C' ; K_{zx} Q_X y = (C^{-1} K_{xz})' (C^{-1} y)
    cholesky(Kxx, lower=True, overwrite_a=True)
    Ci_y = solve_triangular(Kxx, y.copy('F'), lower=True, overwrite_b=True)
## Compute the border vector and the corner value
    Kxz = np.asfortranarray(kernel(X, Z, **kwargs))
    Kxz = solve_triangular(Kxx, Kxz, lower=True, overwrite_b=True)
## Compute the KzxQx product and diagonalize M
    neg_M = np.einsum("ji,ji->i", Kxz, Kxz) - np.diag(kernel(Z, **kwargs)) - nugget
## Start preparing the A matrix:
    y_z = np.dot(Ci_y.T, Kxz)
    solve_triangular(Kxx, Kxz, lower=True, trans=1, overwrite_b=True)
    solve_triangular(Kxx, Ci_y, lower=True, trans=1, overwrite_b=True)
    # Ci_y = Q_X y; Kxz = Q_X K_XZ
## Invert the Cholseky matrix and get the diagonal of L^{-1}' L^{-1}
    dtrtri(Kxx, overwrite_c=1, lower=1)
    diag_ = np.einsum("ji,ji->i", Kxx, Kxx)
    del Kxx
## Now compute the B vector
    B = np.empty((Kxz.shape[0]+1, Kxz.shape[1]))
    B[:-1] = Kxz
    B[-1] = -1.0
    B /= neg_M
## Now it's A'turn: column-wise multiplication of B
    A = - y_z[:, np.newaxis] * B
    A[..., :-1, :] += Ci_y.T[..., np.newaxis]
    A, B = A.swapaxes(-2, -1), B.T
    if use_loo:
        # Kxz = Q_X K_XZ
        Kxz *= -B[:, :-1].T ## NxK
        # Kxz = (Q_X K_XZ)^2 / M_Z
        Kxz += diag_[:, np.newaxis]
## `diag` contains the reciprocals of the leverages (up to `nugget` scaling)
        A[..., :-1] /= Kxz.T
        B[..., :-1] /= Kxz.T
        A[..., -1] *= -neg_M
        B[..., -1] *= -neg_M
    else:
        A *= nugget
        B *= nugget
    del Kxz
## return A, B, KzxQx, M and LOO-residuals: only M depends on sigma2.
    return A, B, y_z.T, - neg_M[:, np.newaxis], Ci_y
