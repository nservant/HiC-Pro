import cython
import numpy as np
cimport numpy as np


cdef extern from "read.c":
    int get_num_line (char* filename)
    void read_counts(char* filename, int* array)

cdef extern from "write.c":
    void write_counts(char* filename, int* x, int* y, double* counts, int n_lines)
    void write_counts_int(char* filename, int* x, int* y, int* counts,
                          int n_lines)


def get_nline(filename):
    n_lines = get_num_line(filename)
    return n_lines


def read(np.ndarray[int, ndim=1, mode="c"] X not None, filename):
    n_lines = get_num_line(filename)
    read_counts(filename, <int*> X.data)
    return n_lines


def loadtxt(filename):
    """
    Fast loading of a raw interaction counts file

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be of the following format:
            i, j, counts

    Returns
    --------
    X : the interaction counts file
    """
    n_lines = get_num_line(filename)
    # We need twice the amount of line because we don't use the fact that the
    # matrix is symmetric, and we need three cols per row
    X = np.empty((n_lines * 3), dtype=np.int32)
    read(X, filename)
    return X.reshape((n_lines, 3))


def savetxt(filename, np.ndarray[int, ndim=1, mode="c"] col not None,
            np.ndarray[int, ndim=1, mode="c"] row not None,
            np.ndarray counts not None):
    """
    Fast writing of interaction counts

    Parameters
    """
    if counts.dtype == float:
        save_txt_float(filename, col, row, counts)
    else:
        save_txt_int(filename, col, row, counts)
    # FIXME


def save_txt_float(filename,
            np.ndarray[int, ndim=1, mode="c"] col not None,
            np.ndarray[int, ndim=1, mode="c"] row not None,
            np.ndarray[double, ndim=1, mode="c"] counts not None):
    """
    Fast writing of normalized interaction counts

    Parameters
    """
    n_lines = col.shape[0]
    write_counts(filename, <int*> row.data, <int*> col.data, <double*> counts.data,
          n_lines)


def save_txt_int(filename,
            np.ndarray[int, ndim=1, mode="c"] col not None,
            np.ndarray[int, ndim=1, mode="c"] row not None,
            np.ndarray[int, ndim=1, mode="c"] counts not None):
    """
    Fast writing of normalized interaction counts

    Parameters
    """
    n_lines = col.shape[0]
    write_counts_int(filename, <int*> row.data, <int*> col.data,
                     <int*> counts.data,
                     n_lines)
