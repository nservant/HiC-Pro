import numpy as np
from nose.tools import assert_raises
from malaria.utils import validation


def test_is_symetric_or_tri():
    n = 100
    m = 50
    X = np.random.random((n, m))
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = np.random.random((n, n))
    assert_raises(ValueError, validation.is_symetric_or_tri, X)
    X = X + X.T
    validation.is_symetric_or_tri(X)
    X[np.tri(n, dtype=bool)] = 0
    validation.is_symetric_or_tri(X)
