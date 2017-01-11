import numpy as np
from scipy import sparse
from .fastio_ import loadtxt


def load_counts(filename, lengths=None):
    """
    Fast loading of a raw interaction counts file

    Parameters
    ----------
    filename : str
        path to the file to load. The file should be of the following format:
        i, j, counts

    lengths : ndarray
        lengths of each chromosomes

    Returns
    --------
    X : the interaction counts file
    """
    n = None
    if lengths is not None:
        n = lengths.sum()
        shape = (n, n)
    else:
        shape = None
    X = loadtxt(filename.encode())
    X = X[X[:, 2] != 0]
    if shape is not None:
        if X[:, :2].max() == shape[0]:
            X[:, :2] -= 1
    counts = sparse.coo_matrix((X[:, 2], (X[:, 0], X[:, 1])),
                               shape=shape, dtype=np.float64)
    return counts


def load_lengths(filename):
    """
    Fast loading of the bed files

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    -------
    lengths : the lengths of each chromosomes
    """
    data = np.loadtxt(filename, dtype="str")
    lengths = [(data[:, 0] == i).sum() for i in np.unique(data[:, 0])]
    return np.array(lengths)


def write_counts(filename, counts):
    """
    Write counts

    Parameters
    ----------

    filename : str

    counts: array-like
    """
    raise NotImplementedError
