import numpy as np
from scipy import sparse
from ._normalization_ import _update_normalization_csr
from .utils import is_symetric_or_tri, is_tri


def ICE_normalization(X, SS=None, max_iter=3000, eps=1e-4, copy=True,
                      norm='l1', verbose=0, output_bias=False,
                      total_counts=None, counts_profile=None):
    """
    ICE normalization

    The imakaev normalization of Hi-C data consists of iteratively estimating
    the bias such that all the rows and columns (ie loci) have equal
    visibility.

    Parameters
    ----------
    X : ndarray or sparse array (n, n)
        raw interaction frequency matrix

    max_iter : integer, optional, default: 3000
        Maximum number of iteration

    eps : float, optional, default: 1e-4
        the relative increment in the results before declaring convergence.

    copy : boolean, optional, default: True
        If copy is True, the original data is not modified.

    norm : string, optional, default: l1
        If set to "l1", will compute the ICE algorithm of the paper. Else, the
        algorithm is adapted to use the l2 norm, as suggested in the SCN
        paper.

    output_bias : boolean, optional, default: False
        whether to output the bias vector.

    total_counts : float, optional, default: None
        the total number of contact counts that the normalized matrix should
        contain. If set to None, the normalized contact count matrix will be
        such that the total number of contact counts equals the initial number
        of interactions.

    Returns
    -------
    X, (bias) : ndarray (n, n)
        Normalized IF matrix and bias of output_bias is True

    Example
    -------
    .. plot:: examples/normalization/plot_ICE_normalization.py
    """
    if copy:
        X = X.copy()

    if sparse.issparse(X):
        if not sparse.isspmatrix_csr(X):
            X = sparse.csr_matrix(X, dtype="float")
        X.sort_indices()
    else:
        X[np.isnan(X)] = 0
    X = X.astype('float')

    m = X.shape[0]
    is_symetric_or_tri(X)
    old_bias = None
    bias = np.ones((m, 1))
    _is_tri = is_tri(X)
    if verbose:
        print("Matrix is triangular superior")

    if total_counts is None:
        total_counts = X.sum()
    for it in np.arange(max_iter):
        if norm == 'l1':
            # Actually, this should be done if the matrix is diag sup or diag
            # inf
            if _is_tri:
                sum_ds = X.sum(axis=0) + X.sum(axis=1).T - X.diagonal()
            else:
                sum_ds = X.sum(axis=0)
        elif norm == 'l2':
            if _is_tri:
                sum_ds = ((X**2).sum(axis=0) +
                          (X**2).sum(axis=1).T -
                          (X**2).diagonal())
            else:
                sum_ds = (X**2).sum(axis=0)

        if SS is not None:
            raise NotImplementedError
        dbias = sum_ds.reshape((m, 1))
        if counts_profile is not None:
            dbias /= counts_profile[:, np.newaxis]
        # To avoid numerical instabilities
        dbias /= dbias[dbias != 0].mean()

        dbias[dbias == 0] = 1
        bias *= dbias

        if sparse.issparse(X):
            X = _update_normalization_csr(X, np.array(dbias).flatten())
        else:
            X /= dbias
            X /= dbias.T

        bias *= np.sqrt(X.sum() / total_counts)
        X *= total_counts / X.sum()

        if old_bias is not None and np.abs(old_bias - bias).sum() < eps:
            if verbose > 1:
                print("break at iteration %d" % (it,))
            break

        if verbose > 1 and old_bias is not None:
            print('ICE at iteration %d %s' %
                  (it, np.abs(old_bias - bias).sum()))
        old_bias = bias.copy()
    # Now that we are finished with the bias estimation, set all biases
    # corresponding to filtered rows to np.nan
    if sparse.issparse(X):
        X = X.tocoo()
        to_rm = (np.array(X.sum(axis=0)).flatten() +
                 np.array(X.sum(axis=1)).flatten()) == 0
    else:
        to_rm = (X.sum(axis=0) + X.sum(axis=1)) == 0
    bias[to_rm] = np.nan
    if output_bias:
        return X, bias
    else:
        return X


def SCN_normalization(X, max_iter=300, eps=1e-6, copy=True):
    """
    Sequential Component Normalization

    Parameters
    ----------
    X : ndarray (n, n)
        raw interaction frequency matrix

    max_iter : integer, optional, default: 300
        Maximum number of iteration

    eps : float, optional, default: 1e-6
        the relative increment in the results before declaring convergence.

    copy : boolean, optional, default: True
        If copy is True, the original data is not modified.

    Returns
    -------
    X : ndarray,
        Normalized IF
    """
    # X needs to be square, else it's gonna fail

    m, n = X.shape
    if m != n:
        raise ValueError

    if copy:
        X = X.copy()
    X = X.astype(float)

    for it in np.arange(max_iter):
        sum_X = np.sqrt((X ** 2).sum(0))
        sum_X[sum_X == 0] = 1
        X /= sum_X
        X = X.T
        sum_X = np.sqrt((X ** 2).sum(0))
        sum_X[sum_X == 0] = 1
        X /= sum_X
        X = X.T

        if np.abs(X - X.T).sum() < eps:
            print("break at iteration %d" % (it,))
            break

    return X
