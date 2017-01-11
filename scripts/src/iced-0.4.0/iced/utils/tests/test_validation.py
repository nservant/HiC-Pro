import numpy as np
from scipy import sparse
from nose.tools import assert_raises
from iced.utils import validation


def test_is_symetric_or_tri():
    n = 100
    m = 50
    random_state = np.random.RandomState(seed=42)
    X = random_state.randn(n, m)
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = random_state.randn(n, n)
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = X + X.T
    validation.is_symetric_or_tri(X)
    X = np.triu(X)
    validation.is_symetric_or_tri(X)


def test_is_symetric_or_tri_sparse():
    n = 100
    m = 50
    random_state = np.random.RandomState(seed=42)
    X = sparse.csr_matrix(random_state.randn(n, m))
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = sparse.csr_matrix(random_state.randn(n, n))
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = random_state.randn(n, n)
    X = X + X.T
    X = sparse.csr_matrix(X)
    validation.is_symetric_or_tri(X)
    X[np.tri(n, dtype=bool)] = 0
    validation.is_symetric_or_tri(X)


def test_is_tri():
    n = 100
    random_state = np.random.RandomState(seed=42)
    X = random_state.randn(n, n)
    assert validation.is_tri(np.triu(X))
    assert validation.is_tri(np.tril(X))
