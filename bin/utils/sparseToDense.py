#! /usr/bin/env python

import sys
import numpy as np
import os
from iced import io
from scipy import sparse

def load_bed(filename):
    '''
    Load a BED file using numpy

    Parameters
    ----------
    filename : str,
        path to the file to load. The file should be a bed file

    Returns
    ------
    data : the tree first columns of the BED file (chromosome, start, end)

    '''
    data = np.genfromtxt(filename, dtype='str')
    return data[:,:3]

def load_lengths_perchr(filename, add_name=True):
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
    u, idx = np.unique(data[:, 0], return_index=True)
    lengths = [(data[:, 0] == i).sum() for i in u[np.argsort(idx)]]
    if add_name:
        return (np.array(lengths), u[np.argsort(idx)])
    else:
        return np.array(lengths)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-b", "--bins", help="BED file with bins coordinates. If provided the chromosome lengths are used to define the output matrix size")
    parser.add_argument("-d", "--di", help="If specified the output matrix is formatted for Dixon et al. TADs calling. In this case --bins is required", action='store_true')
    parser.add_argument("-c", "--perchr", help="If specified intrachromosomal maps are written per chromosome as individual dense matrices. In this case, --bins must also be specified", action='store_true')
    parser.add_argument("-o", "--output", help="Output filename")

    args = parser.parse_args()

    if args.di is True and args.bins is None:
        print "--bins parameter is required when --di is specified"
        sys.exit(1)
        
    if args.perchr is True and args.bins is None:
        print "--bins parameter is required when --perchr is specified"
        sys.exit(1)
    
    ## bin option
    if args.bins is not None:
        chr_lengths = load_lengths_perchr(args.bins)
        lengths = chr_lengths[0]
        chrnames = chr_lengths[1]
    else:
        lengths = None

    ## Load counts in sparse format
    counts = io.load_counts(args.filename, lengths=lengths)
 
    ## di option
    if args.di is True:
        bins = load_bed(args.bins)
        if len(bins) != counts.shape[1]:
            print "Error -  number of rows in BED and matrix files are not equal"
            sys.exit(1)
    
    ## Genome-wide dense matrix
    if args.perchr is False:
        counts = counts.toarray()
        counts = counts + counts.T
        counts[np.diag_indices_from(counts)] /= 2
        counts = np.round(counts, 3)

        ## Output name for save
        if args.output is None:
            output_name = os.path.basename(args.filename)
            output_name = output_name.replace(".matrix", "_dense.matrix")
        else:
            output_name = args.output
        
        if args.di is True:
            counts = np.hstack((bins, counts))

        ## Save
        np.savetxt(output_name, counts, '%s', delimiter="\t")
        
    ## Per chromosome matrices
    else:
        ## indexes of intra chrom maps
        lc = np.concatenate([np.array([0]), lengths.cumsum()])
        #print lc
        #print lengths

        for i in range(1, len(lc)):
            print str(chrnames[i-1]) + "..."
            idxintra = np.where(((counts.row >= lc[i-1]) & (counts.row<lc[i])) & ((counts.col>=lc[i-1]) & (counts.col<lc[i])))[0]
         
            ## Subset the counts array and rescale the index based on cumulative lengths
            counts_perchr = sparse.coo_matrix((counts.data[idxintra], (counts.row[idxintra] - lc[i-1], counts.col[idxintra] - lc[i-1])), shape=(lengths[i-1], lengths[i-1]))
            counts_perchr = counts_perchr.toarray()
            counts_perchr = counts_perchr + counts_perchr.T
            counts_perchr[np.diag_indices_from(counts_perchr)] /= 2
            counts_perchr = np.round(counts_perchr, 3)
            
            ## Output name for save
            if args.output is None:
                output_name = args.filename.replace(".matrix", "_" + str(chrnames[i-1]) + "_dense.matrix")
            else:
                output_name = str(chrnames[i-1]) + "_" + args.output 
        
            if args.di is True:
                counts_perchr = np.hstack((bins[np.where(bins==str(chrnames[i-1]))[0]], counts_perchr))

            ## Save
            np.savetxt(output_name, counts_perchr, '%s', delimiter="\t")
