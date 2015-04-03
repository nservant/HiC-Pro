import numpy as np
from scipy import sparse


def is_symetric_or_tri(X, eps=1e-7):
    m, n = X.shape
    if m != n:
        raise ValueError("The matrix should be of shape (n, n)")
    if sparse.issparse(X):
        return True
    if (X - X.T >= eps).any() and X[np.tri(n, dtype=bool)].any():
        raise ValueError("The matrix should be symmetric")
