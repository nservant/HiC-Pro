from iced.filter import filter_low_counts, filter_high_counts
import numpy as np
from numpy.testing import assert_array_equal
from scipy import sparse


def test_filter_low_counts():
    X = np.ones((100, 100))
    X[0, :] = 0
    X[:, 0] = 0
    X_filtered_true = X.copy()
    X_filtered_true[X == 0] = np.nan
    X_filtered = filter_low_counts(X)
    assert_array_equal(X_filtered, X_filtered_true)

    lengths = np.array([40, 60])
    X_filtered = filter_low_counts(X, lengths=lengths)
    assert_array_equal(X_filtered, X_filtered_true)

    X_filtered = filter_low_counts(X, sparsity=False)
    assert_array_equal(X_filtered, X_filtered_true)


def test_sparse_filter_low_counts():
    X = 10 * np.ones((100, 100))
    X[0, :] = 1
    X[:, 0] = 1
    X_filtered_dense = X.copy()
    X_filtered_dense[0] = 0
    X_filtered_dense[:, 0] = 0
    return
    # this is not implemented yet
    X_filtered_sparse_csr = filter_low_counts(sparse.csr_matrix(X),
                                              sparsity=False)
    X_filtered_sparse_coo = filter_low_counts(sparse.coo_matrix(X))

    assert_array_equal(X_filtered_dense,
                       np.array(X_filtered_sparse_csr.todense()))
    assert_array_equal(X_filtered_dense,
                       np.array(X_filtered_sparse_coo.todense()))


def test_filter_high_counts():
    X = np.ones((100, 100))
    X *= np.arange(len(X))
    X_filtered_true = X.copy()
    X_filtered_true[-1] = np.nan
    X_filtered_true[:, -1] = np.nan
    X_filtered = filter_high_counts(X)
    assert_array_equal(X_filtered, X_filtered_true)


def test_filter_high_couts_sparse():
    X = np.ones((100, 100))
    X *= np.arange(len(X))
    X_filtered_dense = filter_high_counts(X)

    X_filtered_sparse_csr = filter_high_counts(sparse.csr_matrix(X))
    X_filtered_sparse_coo = filter_high_counts(sparse.coo_matrix(X))

    assert_array_equal(X_filtered_dense,
                       np.array(X_filtered_sparse_csr.todense()))
    assert_array_equal(X_filtered_dense,
                       np.array(X_filtered_sparse_coo.todense()))
