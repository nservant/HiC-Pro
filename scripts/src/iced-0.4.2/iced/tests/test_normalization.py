import numpy as np
from scipy import sparse
from numpy.testing import assert_array_almost_equal
import nose
from nose.tools import assert_almost_equal

from iced.normalization import ICE_normalization
from iced.normalization import SCN_normalization


def test_ICE_normalization():
    n = 100
    random_state = np.random.RandomState(seed=42)
    X = random_state.randint(0, 100, size=(n, n))
    X = X + X.T
    normed_X = ICE_normalization(X, eps=1e-10, max_iter=1000000)
    normed = normed_X.sum(axis=1)
    assert_array_almost_equal(normed / normed.mean(), np.ones((len(X), )),
                              decimal=0)

    normed_X, bias = ICE_normalization(X, eps=1e-10, max_iter=100000,
                                       output_bias=True)
    assert_array_almost_equal(normed_X, X / (bias.T * bias), 6)


def test_sparse_ICE_normalization():
    n = 100
    random_state = np.random.RandomState(seed=42)
    X = random_state.randint(0, 100, size=(n, n))

    thres = (random_state.rand(n, n) > 0.5).astype(bool)

    X[thres] = 0
    X = X + X.T
    sparse_X = sparse.csr_matrix(X)
    true_normed_X = ICE_normalization(X, eps=1e-10, max_iter=10)
    normed_X = ICE_normalization(sparse_X, eps=1e-10, max_iter=10)
    assert_array_almost_equal(X, sparse_X.todense())
    assert_array_almost_equal(true_normed_X, np.array(normed_X.todense()))


def test_sparse_ICE_normalization_triu():
    n = 100
    random_state = np.random.RandomState(seed=42)
    X = random_state.randint(0, 100, size=(n, n))

    thres = (random_state.rand(n, n) > 0.5).astype(bool)
    X[thres] = 0
    X = X + X.T
    sparse_X = sparse.triu(X)
    true_normed_X, true_biases = ICE_normalization(
        X, eps=1e-10, max_iter=10, output_bias=True)
    true_normed_X = np.triu(true_normed_X)

    normed_X_sparse, biases_sparse = ICE_normalization(
        sparse_X, eps=1e-10, max_iter=100,
        output_bias=True)
    normed_X_dense, biases_dense = ICE_normalization(
        np.triu(X), eps=1e-10, max_iter=100,
        output_bias=True)

    # The sparse and dense version are going to be equal up to a constant
    # factor
    assert_array_almost_equal(normed_X_dense,
                              np.array(normed_X_sparse.toarray()))

    normed_X_sparse *= true_normed_X.mean() / normed_X_sparse.mean()
    normed_X_dense *= true_normed_X.mean() / normed_X_dense.mean()

    assert_array_almost_equal(true_normed_X,
                              np.array(normed_X_sparse.todense()))
    assert_array_almost_equal(true_normed_X, normed_X_dense)

    total_counts = 5000
    normed_X = ICE_normalization(sparse_X, eps=1e-10,
                                 total_counts=total_counts)
    assert_almost_equal(normed_X.sum(), total_counts)


def test_SCN_normalization():
    n = 100
    random_state = np.random.RandomState(seed=42)
    X = random_state.randint(0, 100, size=(n, n))

    normed_X = SCN_normalization(X)
    assert_array_almost_equal(np.sqrt((normed_X ** 2).sum(axis=1)),
                              np.ones((len(X), )))


if __name__ == "__main__":
    nose.runmodule(argv=['-s', '--with-doctest'], exit=False)
