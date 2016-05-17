#! /usr/bin/env python

import numpy as np
from iced import io


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--lengths", "-l")
    parser.add_argument("--output", "-o")

    args = parser.parse_args()

    if args.lengths is not None:
        lengths = io.load_lengths(args.lengths)
    else:
        lengths = None

    counts = io.load_counts(args.filename, lengths=lengths)
    counts = counts.toarray()
    counts = counts + counts.T
    counts[np.diag_indices_from(counts)] /= 2

    # save matrix like file
    if args.output is None:
        output_name = args.filename.replace(".matrix", "_dense.matrix")
    else:
        output_name = args.output

    np.savetxt(output_name, counts, '%s')
