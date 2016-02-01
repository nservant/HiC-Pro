import numpy as np
cimport cython
cimport numpy as np

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INT


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _filter_csr(X, np.ndarray[INT, ndim=1] bias):

    cdef:
        np.ndarray[DOUBLE, ndim=1] X_data = X.data
        np.ndarray[int, ndim=1] X_indices = X.indices
        np.ndarray[int, ndim=1] X_indptr = X.indptr
        unsigned int m = X.shape[0]
        unsigned int i, j, row

    j = 0
    for i, row in enumerate(X_indices):
        if i >= X_indptr[j + 1]:
            j += 1
        if bias[row] or bias[j]:
            X_data[i] = 0
    return X
