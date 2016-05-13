#! /usr/bin/env python

import sys
import numpy as np
from scipy import sparse


def open_with_numpy_loadtxt(filename):
    '''
    http://stackoverflow.com/questions/4315506/load-csv-into-2d-matrix-with-numpy-for-plotting
    '''
    data = np.loadtxt(open(filename,'rb'),delimiter="\t",skiprows=0)
    return data

def open_with_pandas_read_csv(filename):
    import pandas as pd
    df = pd.read_csv(filename, sep="\t", header=None)
    data = df.values
    return data    

def open_with(filename):
    ## replace by try / catch import pandas
    return open_with_numpy_loadtxt(filename)

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
    # This is the interaction count files
    df = open_with(filename)
    row, col, data = df.T
    # XXX We need to deal with the fact that we should not duplicate entries
    # for the diagonal.
    # XXX what if n doesn't exist?
    if (col.min() >= 1 and row.min() >= 1) and \
       ((n is None) or (col.max() == n)):
        # This is a hack to deal with the fact that sometimes, the files are
        # indexed at 1 and not 0
        col -= 1
        row -= 1

    data = data.astype(float)
    counts = sparse.coo_matrix((data, (row, col)), shape=shape)
    return counts


def load_lengths(filename):
    """
    Fast loading of the chromosome length files

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    -------
    lengths : the lengths of each chromosomes
    """
    data = open_with_pandas_read_csv(filename)
    lengths = [(data[:, 0] == i).sum() for i in np.unique(data[:, 0])]
    return np.array(lengths)

                                                                                                                                                                         
def load_bed(filename):                                                                                                                
    data = open_with_pandas_read_csv(filename)
    return data[:,:3]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--lengths", "-l")
    parser.add_argument("--bins", "-b")
    parser.add_argument("--output", "-o")
    
    args = parser.parse_args()

    if args.lengths is not None:
        lengths = load_lengths(args.lengths)
    else:
        lengths = None

    print lengths

    counts = load_counts(args.filename, lengths=lengths)
    counts = counts.toarray()
    counts = counts + counts.T
    counts[np.diag_indices_from(counts)] /= 2
    
    if args.bins is not None:
        bins = load_bed(args.bins)
        if len(bins) != len(counts):
            print "Error -  number of rows in BED and matrix files are not equal"
            sys.exit(1)
        counts = np.hstack((bins, counts))
    else:
        bins = None

    ## save matrix like file
    if args.output is None:
        outputName = args.filename.replace(".matrix","_dense.matrix")
    else:
        outputName = args.output
    
    np.savetxt(outputName, counts, '%s')
  
