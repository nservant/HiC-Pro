import os
from os import environ, makedirs
from os.path import join, expanduser, exists
from os.path import dirname
import shutil

from .. import io

# authors: Nelle Varoquaux <nelle.varoquaux@gmail.com>

# This module is greatly inspired from sklearn.datasets


def get_data_home(data_home=None):
    """Return the path of the scikit-learn data dir.

    This folder is used by some large dataset loaders to avoid
    downloading the data several times.

    By default the data dir is set to a folder named 'scikit_learn_data'
    in the user home folder.

    Alternatively, it can be set by the 'SCIKIT_LEARN_DATA' environment
    variable or programmatically by giving an explicit folder path. The
    '~' symbol is expanded to the user home folder.

    If the folder does not already exist, it is automatically created.
    """
    if data_home is None:
        data_home = environ.get('HICLIB_DATA',
                                join('~', 'hiclib_data'))
    data_home = expanduser(data_home)
    if not exists(data_home):
        makedirs(data_home)
    return data_home


def clear_data_home(data_home=None):
    """Delete all the content of the data home cache."""
    data_home = get_data_home(data_home)
    shutil.rmtree(data_home)


def load_sample_yeast():
    """
    Load and return a sample of S. cerevisiae contact count matrix from duan
    et al, Nature, 2009

    Returns
    -------
        counts, lengths:
            tuple of two elements, the first a contact count matrix, the
            second an ndarray containing the lengths of the chromosomes.
    """
    module_path = dirname(__file__)
    lengths = io.load_lengths(
        os.path.join(module_path, "data/duan2009/duan.SC.10000.raw_sub.bed"))
    counts = io.load_counts(
        os.path.join(module_path,
                     "data/duan2009/duan.SC.10000.raw_sub.matrix"),
        lengths=lengths)
    counts = counts.toarray()
    counts = counts.T + counts
    return counts, lengths
