#! /usr/bin/env python

import sys
import argparse
import numpy as np
from scipy import sparse

import iced
from iced.io import loadtxt, savetxt


parser = argparse.ArgumentParser("ICE normalization")
parser.add_argument('filename',
                    metavar='File to load',
                    type=str,
                    help='Path to file of contact counts to load')
parser.add_argument("--results_filename",
                    "-r",
                    type=str,
                    default=None,
                    help="results_filename")
parser.add_argument("--filtering_perc", "-f",
                    type=float,
                    default=None,
                    help="Percentage of reads to filter out")
parser.add_argument("--filter_low_counts_perc",
                    type=float,
                    default=0.02,
                    help="Percentage of reads to filter out")
parser.add_argument("--filter_high_counts_perc",
                    type=float,
                    default=0,
                    help="Percentage of reads to filter out")
parser.add_argument("--remove-all-zeros-loci", default=False,
                    action="store_true",
                    help="If provided, all non-interacting loci will be "
                         "removed prior to the filtering strategy chosen.")
parser.add_argument("--max_iter", "-m", default=100, type=int,
                    help="Maximum number of iterations")
parser.add_argument("--eps", "-e", default=0.1, type=float,
                    help="Precision")
parser.add_argument("--dense", "-d", default=False, action="store_true")
parser.add_argument("--output-bias", "-b", default=False, help="Output the bias vector")
parser.add_argument("--verbose", "-v", default=False)


args = parser.parse_args()
filename = args.filename

# Deprecating filtering_perc option
filter_low_counts = None
if "--filtering_perc" in sys.argv:
    DeprecationWarning(
        "Option '--filtering_perc' is deprecated. Please use "
        "'--filter_low_counts_perc' instead.'")
    # And print it again because deprecation warnings are not displayed for
    # recent versions of python
    print "--filtering_perc is deprecated. Please use filter_low_counts_perc"
    print "instead. This option will be removed in ice 0.3"
    filter_low_counts = args.filtering_perc
if "--filter_low_counts_perc" in sys.argv and "--filtering_perc" in sys.argv:
    raise Warning("This two options are incompatible")
if "--filtering_perc" is None and "--filter_low_counts_perc" not in sys.argv:
    filter_low_counts_perc = 0.02
elif args.filter_low_counts_perc is not None:
    filter_low_counts_perc = args.filter_low_counts_perc

if args.verbose:
    print("Using iced version %s" % iced.__version__)
    print "Loading files..."

# Loads file as i, j, counts
i, j, data = loadtxt(filename).T

# Detecting whether the file is 0 or 1 based.
if min(i.min(), j.min()) == 0:
    index_base = 0
    N = max(i.max(), j.max()) + 1
    counts = sparse.coo_matrix((data, (i, j)), shape=(N, N), dtype=float)
else:
    index_base = 1
    N = max(i.max(), j.max())
    counts = sparse.coo_matrix((data, (i - 1, j - 1)), shape=(N, N), dtype=float)

if args.dense:
    counts = np.array(counts.todense())
else:
    counts = sparse.csr_matrix(counts)

if args.verbose:
    print "Normalizing..."

if filter_low_counts_perc != 0:
    counts = iced.filter.filter_low_counts(counts,
                                           percentage=filter_low_counts_perc,
                                           remove_all_zeros_loci=args.remove_all_zeros_loci,
                                           copy=False, sparsity=False, verbose=args.verbose)
if args.filter_high_counts_perc != 0:
    counts = iced.filter.filter_high_counts(
        counts,
        percentage=args.filter_high_counts_perc,
        copy=False)

counts, bias = iced.normalization.ICE_normalization(
    counts, max_iter=args.max_iter, copy=False,
    verbose=args.verbose, eps=args.eps, output_bias=True)

if args.results_filename is None:
    results_filename = ".".join(
        filename.split(".")[:-1]) + "_normalized." + filename.split(".")[-1]
else:
    results_filename = args.results_filename

counts = sparse.coo_matrix(counts)

if args.verbose:
    print "Writing results..."

savetxt(
    results_filename, counts.col + index_base, counts.row + index_base, counts.data)


if args.output_bias:
    np.savetxt(results_filename + ".biases", bias)
