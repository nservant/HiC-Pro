import numpy as np
from scipy import sparse


def is_symetric_or_tri(X, eps=1e-7):
    m, n = X.shape
    if m != n:
        raise ValueError("The matrix should be of shape (n, n)")

    if is_tri(X):
        return True
    if np.abs(X - X.T).sum() > eps:
        raise ValueError("The matrix should be symmetric")


def is_tri(X):
    diag = X.diagonal().sum()
    if sparse.issparse(X):
        if not (sparse.tril(X).sum() - diag) or \
           not (sparse.triu(X).sum() - diag):
            return True
    elif not np.triu(X, 1).sum() or not np.tril(X, -1).sum():
        return True
    else:
        return False
