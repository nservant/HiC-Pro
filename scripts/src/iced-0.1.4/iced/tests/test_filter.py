from iced._filter import filter_low_counts
import numpy as np
from numpy.testing import assert_array_equal
from scipy import sparse


def test_filter_low_counts():
    X = np.ones((100, 100))
    X[0, :] = 0
    X[:, 0] = 0
    X_filtered = filter_low_counts(X)
    X_filtered_true = X.copy()
    X_filtered_true[X == 0] = np.nan
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
    X_filtered = filter_low_counts(sparse.csr_matrix(X), sparsity=False)
    X_filtered_true = X.copy()
    X_filtered_true[0, :] = 0
    X_filtered_true[:, 0] = 0
    assert_array_equal(X_filtered.todense(),
                       X_filtered_true)
