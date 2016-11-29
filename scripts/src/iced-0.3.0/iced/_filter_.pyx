import numpy as np
cimport cython
cimport numpy as cnp

ctypedef cnp.float64_t DOUBLE
ctypedef cnp.int8_t BOOL


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _filter_csr(X, cnp.ndarray[BOOL, ndim=1, cast=True] bias):

    cdef:
        cnp.ndarray[DOUBLE, ndim=1] X_data = X.data
        cnp.ndarray[int, ndim=1] X_indices = X.indices
        cnp.ndarray[int, ndim=1] X_indptr = X.indptr
        unsigned int m = X.shape[0]
        unsigned int i, j, row

    j = 0
    for i, row in enumerate(X_indices):
        while i >= X_indptr[j + 1]:
            j += 1
        if bias[row] or bias[j]:
            X_data[i] = 0
    return X
