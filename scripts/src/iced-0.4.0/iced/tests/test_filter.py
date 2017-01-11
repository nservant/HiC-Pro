from iced.filter import filter_low_counts, filter_high_counts
from iced.datasets import load_sample_yeast
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


def test_filter_low_counts_with_weros():
    X = 10 * np.ones((100, 100))
    X[0, :] = 0
    X[:, 0] = 0
    X[1, :] = 1
    X[:, 1] = 1

    X_filtered_true = X.copy()
    X_filtered_true[X != 10] = np.nan
    X_filtered = filter_low_counts(X, remove_all_zeros_loci=True,
                                   sparsity=False)
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


def test_sparse_filter_low_counts_real_data():
    counts, lengths = load_sample_yeast()
    counts_sparse = sparse.csr_matrix(counts)
    counts_dense = filter_low_counts(counts, sparsity=False, percentage=0.1)
    counts_sparse = filter_low_counts(counts_sparse, sparsity=False,
                                      percentage=0.1)
    counts_dense[np.isnan(counts_dense)] = 0
    assert_array_equal(counts_dense, counts_sparse.toarray())

    triu_counts_sparse = sparse.csr_matrix(np.triu(counts))
    triu_counts_sparse = filter_low_counts(triu_counts_sparse, sparsity=False,
                                           percentage=0.1)
    assert_array_equal(np.triu(counts), triu_counts_sparse.toarray())


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
