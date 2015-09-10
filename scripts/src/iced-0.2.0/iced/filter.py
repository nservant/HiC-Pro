import numpy as np
from scipy import sparse
from ._filter_ import _filter_csr
from . import utils


def filter_low_counts(X, lengths=None, percentage=0.02, copy=True,
                      sparsity=True):
    """
    Filter rows and columns with low counts

    Parameters
    ----------
    X : ndarray (n, n)
        Count matrix (hollow, symetric)

    lengths : ndarray (m, ), optional, default: None
        Lengths of the chromosomes

    percentage : float, optional, default: 0.02
        percentage of rows and columns to discard

    copy : boolean, optional, default: True
        If set to true, copies the count matrix

    sparsity : boolean, optional, default: True
        Whether to use the sparsity of the region or the total number of
        counts of the region to filer.

    Return
    ------
    X : ndarray (n, n)
        The filtered array
    """
    if sparse.issparse(X):
        if not sparse.isspmatrix_csr(X):
            X = X.tocsr()
        X.sort_indices()
    else:
        X[np.isnan(X)] = 0

    if sparsity:
        if lengths is not None:
            weights = []
            [weights.append(i) for i in lengths for j in range(i)]
            weights = np.array(weights)
            mask = utils.get_inter_mask(lengths, resolution=1)
        else:
            weights = np.ones(X.shape[0])
            mask = np.zeros(X.shape, dtype=np.bool)

        return _filter_low_sparse(X, weights, mask, percentage=percentage)
    else:
        return _filter_low_sum(X, percentage=percentage)


def filter_high_counts(X, lengths=None, percentage=0.02, copy=True):
    """
    Filter rows and columns with high counts

    Parameters
    ----------
    X : ndarray (n, n)
        Count matrix (hollow, symetric)

    lengths : ndarray (m, ), optional, default: None
        Lengths of the chromosomes

    percentage : float, optional, default: 0.02
        percentage of rows and columns to discard

    copy : boolean, optional, default: True
        If set to true, copies the count matrix

    Return
    ------
    X : ndarray (n, n)
        The filtered array

    Notes
    -----
    New in 0.6
    """
    if sparse.issparse(X):
        if not sparse.isspmatrix_csr(X):
            X = X.tocsr()
        X.sort_indices()
    else:
        X[np.isnan(X)] = 0

    return _filter_high_sum(X, percentage=percentage)


def _filter_low_sparse(X, weights, mask, percentage=0.02):
    # This is NOT going to work on sparse data. For now, raise a Not
    # implemented error

    if sparse.issparse(X):
        raise NotImplemented
    counts = X.copy()
    counts[mask] = 1
    X_sum = (counts == 0).sum(axis=0).astype(float) / weights
    X_sum.sort()
    X_sum = np.array(X_sum).flatten()
    x = X_sum[int(X.shape[0] * (1. - percentage))]
    X_sum = (counts == 0).sum(axis=0).astype(float) / weights

    if sparse.issparse(X):
        _filter_csr(X, (X_sum < x))
    else:
        X[X_sum > x, :] = np.nan
        X[:, X_sum > x] = np.nan

    return X


def _filter_high_sum(X, percentage=0.02):
    X_sum = np.array(X.sum(axis=0)).flatten()
    X_sum.sort()
    m = X.shape[0]
    x = X_sum[int(m * (1-percentage))]
    X_sum = np.array(X.sum(axis=0)).flatten()

    if sparse.issparse(X):
        _filter_csr(X, (X_sum > x).astype(np.int32))
    else:
        X[X_sum > x, :] = np.nan
        X[:, X_sum > x] = np.nan

    return X


def _filter_low_sum(X, percentage=0.02):
    X_sum = np.array(X.sum(axis=0)).flatten()
    X_sum.sort()
    m = X.shape[0]
    x = X_sum[int(m * percentage)]
    X_sum = np.array(X.sum(axis=0)).flatten()

    if sparse.issparse(X):
        _filter_csr(X, (X_sum < x).astype(np.int32))
    else:
        X[X_sum < x, :] = np.nan
        X[:, X_sum < x] = np.nan

    return X
