#! /usr/bin/env python

import sys
import numpy as np
from iced import io

def load_bed(filename):
    '''
    Load a BED file using numpy
    '''
    data = np.genfromtxt(filename, dtype='str')
    return data[:,:3]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-b", "--bins", help="BED file with bins coordinates. If provided the chromosome lengths are used to define the output matrix size")
    parser.add_argument("-d", "--di", help="If specified the output matrix is formatted for Dixon et al. TADs calling. In this case --bins is required", action='store_true')
    parser.add_argument("-o", "--output", help="Output filename")

    args = parser.parse_args()

    if args.di is True and args.bins is None:
        print "--bins parameter is required when --di is specified"
        sys.exit(1)
        
    if args.bins is not None:
        lengths = io.load_lengths(args.bins)
    else:
        lengths = None

    counts = io.load_counts(args.filename, lengths=lengths)
    counts = counts.toarray()
    counts = counts + counts.T
    counts[np.diag_indices_from(counts)] /= 2

    ## rounds matrix
    counts = np.round(counts, 3)

    if args.di is True:
        bins = load_bed(args.bins)
        if len(bins) != len(counts):
            print "Error -  number of rows in BED and matrix files are not equal"
            sys.exit(1)
        counts = np.hstack((bins, counts))

    # save matrix like file
    if args.output is None:
        output_name = args.filename.replace(".matrix", "_dense.matrix")
    else:
        output_name = args.output

    np.savetxt(output_name, counts, '%s', delimiter="\t")
