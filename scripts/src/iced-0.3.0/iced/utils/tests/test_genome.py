import numpy as np
from nose.tools import assert_equal
from numpy.testing import assert_array_equal

from iced.utils._genome import get_intra_mask
from iced.utils._genome import get_inter_mask
from iced.utils._genome import _change_lengths_resolution
from iced.utils._genome import undersample_per_chr
from iced.utils._genome import extract_sub_contact_map
from iced.utils._genome import downsample_resolution


def test_get_intra_mask():
    lengths = np.array([5, 5])
    mask = get_intra_mask(lengths)
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
    mask = get_inter_mask(lengths)
    true_mask = np.zeros((10, 10))
    true_mask[:5, :5] = 1
    true_mask[5:, 5:] = 1
    assert_array_equal(mask, np.invert(true_mask.astype(bool)))


def test_downsample_resolution():
    random_state = np.random.RandomState(seed=42)

    lengths = np.array([10, 10])
    counts = np.triu(random_state.randint(
        0, 100, (lengths.sum(), lengths.sum())))
    counts = counts + counts.T
    downsampled_counts, downsampled_lengths = downsample_resolution(
        counts, lengths)
    assert_equal(downsampled_lengths.sum(), lengths.sum()/2)

    lengths = np.array([10, 11])
    counts = np.triu(random_state.randint(
        0, 100, (lengths.sum(), lengths.sum())))
    counts = counts + counts.T
    downsampled_counts, downsampled_lengths = downsample_resolution(
        counts, lengths)
    assert_equal(downsampled_lengths.sum(), 11)


def test_undersample_per_chr():
    X = np.array([[1, 1, 0, 0],
                  [1, 1, 0, 0],
                  [0, 0, 0.5, 0.5],
                  [0, 0, 0.5, 0.5]])
    lengths = np.array([2, 2])
    undersampled_X = undersample_per_chr(X, lengths)
    undersampled_X_true = np.array([[1, 0],
                                    [0, 0.5]])
    assert_array_equal(undersampled_X_true, undersampled_X)


def test_return_sample():
    lengths = np.array([50, 75])
    n = lengths.sum()
    X = np.random.randint(0, 50, (n, n))
    X = np.triu(X)
    sub_X, _ = extract_sub_contact_map(X, lengths, [0])
    assert_array_equal(X[:lengths[0], :lengths[0]],
                       sub_X)
