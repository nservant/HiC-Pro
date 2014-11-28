import numpy as np
from numpy.testing import assert_array_equal

from malaria.utils._genome import get_intra_mask
from malaria.utils._genome import get_inter_mask
from malaria.utils._genome import _change_lengths_resolution
from malaria.utils._genome import undersample_per_chr


def test_get_intra_mask():
    lengths = np.array([5, 5])
    mask = get_intra_mask(lengths, resolution=1)
    true_mask = np.zeros((10, 10))
    true_mask[:5, :5] = 1
    true_mask[5:, 5:] = 1
    assert_array_equal(mask, true_mask.astype(bool))


def test_change_lengths_resolution():
    lengths = np.array([5, 5])
    l = _change_lengths_resolution(lengths, resolution=1)
    assert_array_equal(lengths, l)


def test_get_inter_mask():
    lengths = np.array([5, 5])
    mask = get_inter_mask(lengths, resolution=1)
    true_mask = np.zeros((10, 10))
    true_mask[:5, :5] = 1
    true_mask[5:, 5:] = 1
    assert_array_equal(mask, true_mask.astype(bool) == False)


def test_undersample_per_chr():
    X = np.array([[1, 1, 0, 0],
                  [1, 1, 0, 0],
                  [0, 0, 0.5, 0.5],
                  [0, 0, 0.5, 0.5]])
    lengths = np.array([2, 2])
    undersampled_X = undersample_per_chr(X, lengths, resolution=1)
    undersampled_X_true = np.array([[1, 0],
                                    [0, 0.5]])
    assert_array_equal(undersampled_X_true, undersampled_X)
