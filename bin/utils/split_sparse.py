#! /usr/bin/env python

import sys
import re
import numpy as np
import os
import iced
from iced.io import loadtxt, load_counts, savetxt
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
    return np.array(lengths)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Split a sparse genome-wide contact map into per chromosome sparse contact map(s)")
    parser.add_argument("filename")
    parser.add_argument("-b", "--bins", help="BED file with bins coordinates. If provided the chromosome lengths are used to define the output matrix size", required=True)
    parser.add_argument("-c", "--chr", help="Specified an intrachromosomal map to export")
    parser.add_argument("-o", "--output", help="Output prefix")

    args = parser.parse_args()
        
    ## bin option
    if args.bins is not None:
        chr_lengths = load_lengths_perchr(args.bins)
        lengths = chr_lengths[0]
        chrnames = chr_lengths[1]

    ## Load bed
    bed = load_bed(args.bins)
        
    ## Load counts in sparse format
    counts = load_counts(args.filename, lengths=lengths)
 
    ## indexes of intra chrom maps
    lc = np.concatenate([np.array([0]), lengths.cumsum()])
    
    for i in range(1, len(lc)):
        if args.chr is None or (args.chr is not None and str(chrnames[i-1]) == args.chr):
            print(str(chrnames[i-1]) + "...")
            idxintra = np.where(((counts.row >= lc[i-1]) & (counts.row<lc[i])) & ((counts.col>=lc[i-1]) & (counts.col<lc[i])))[0]
         
            ## Subset the counts array and rescale the index based on cumulative lengths
            counts_perchr = sparse.coo_matrix((counts.data[idxintra], (counts.row[idxintra] - lc[i-1], counts.col[idxintra] - lc[i-1])), shape=(lengths[i-1], lengths[i-1]))
            
            ## Output name for save
            if args.output is None:
                output_prefix = re.sub(".mat(rix)*", "_" + str(chrnames[i-1]), os.path.basename(args.filename))
            else:
                output_prefix = args.output + "_" + str(chrnames[i-1]) 
            
            ## Save
            mout=np.column_stack((counts_perchr.row + 1, counts_perchr.col + 1,  np.round(counts_perchr.data, 3)))
            savetxt(output_prefix + ".matrix", mout, fmt='%d\t%d\t%.3f')

            ## Bed file
            bedchr = bed[lc[i-1]:lc[i],:]
            nm = np.array(range(1, bedchr.shape[0] + 1)).reshape(bedchr.shape[0],1)
            bedchr = np.append(bedchr, nm, axis=1)
            np.savetxt(output_prefix + "_abs.bed", bedchr, "%s", delimiter="\t")
            
