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
    return np.array(lengths)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-b", "--bins", help="BED file with bins coordinates.\
                         If provided the chromosome lengths are used to\
                         define the output matrix size")
    parser.add_argument("-g", "--org", help="Reference genome.\
                         Used if --ins is specified", default='org')

    parser.add_argument("-d", "--di", help="If specified the output matrix is\
                         formatted for Dixon et al. TADs calling (directionality index).\
                         In this case --bins is required", action='store_true')
    parser.add_argument("-i", "--ins", help="If specified the output matrix is\
                         formatted for Crane et al. TADs calling (insulation score)\
                         .In this case --bins is required", action='store_true')
    parser.add_argument("-c", "--perchr", help="If specified intrachromosomal\
                         maps are written per chromosome as individual dense\
                         matrices. In this case, --bins must also be specified", 
                        action='store_true')
    parser.add_argument("-o", "--output", help="Output filename")

    args = parser.parse_args()

    if args.di is True and args.bins is None:
        print("--bins parameter is required when --di is specified")
        sys.exit(1)

    if args.ins is True and args.bins is None:
        print("--bins parameter is required when --is is specified")
        sys.exit(1)

    if args.perchr is True and args.bins is None:
        print("--bins parameter is required when --perchr is specified")
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

    ## transform to integer if possible
    if counts.data[0].is_integer():
        counts.data = counts.data.astype(int)
    
    ## di/is option
    if args.di is True or args.ins is True:
        bins = load_bed(args.bins)
        if len(bins) != counts.shape[1]:
            print("Error -  number of rows in BED and matrix files are not equal")
            sys.exit(1)

    if args.ins is True:
        def myfunc( x, idx, org):
            return "bin" + str(idx) + "|" + org + "|" + x[0] + ":" + x[1] + "-" + x[2]
        bins_ins = np.array([myfunc(v, i, args.org) for i,v in enumerate(bins)])
        bins_ins = bins_ins.reshape(len(bins_ins), 1)
        
    ## Genome-wide dense matrix
    if args.perchr is False:
        counts = counts.toarray()
        counts = counts + counts.T
        counts[np.diag_indices_from(counts)] = counts[np.diag_indices_from(counts)] / 2
        counts = np.round(counts, 3)

        ## Output name for save
        if args.output is None:
            output_name = os.path.basename(args.filename)
            output_name = output_name.replace(".matrix", "_dense.matrix")
        else:
            output_name = args.output
        
        if args.di is True:
            counts = np.hstack((bins, counts))
        elif args.ins is True:
            ## add col header
            counts = np.hstack((bins_ins, counts))
            ## add row header
            counts = np.vstack((np.insert(bins_ins, 0, ''), counts))
            
        ## Save
        np.savetxt(output_name, counts, '%s', delimiter="\t")
        
    ## Per chromosome matrices
    else:
        ## indexes of intra chrom maps
        lc = np.concatenate([np.array([0]), lengths.cumsum()])

        for i in range(1, len(lc)):
            print("{}...".format(str(chrnames[i-1])))
            idxintra = np.where(((counts.row >= lc[i-1]) & (counts.row<lc[i])) & 
                                ((counts.col>=lc[i-1]) & (counts.col<lc[i])))[0]
         
            ## Subset the counts array and rescale the index based on cumulative lengths
            counts_perchr = sparse.coo_matrix((counts.data[idxintra],
           (counts.row[idxintra] - lc[i-1], counts.col[idxintra] - lc[i-1])), 
            shape=(lengths[i-1], lengths[i-1]))
            counts_perchr = counts_perchr.toarray()
            counts_perchr = counts_perchr + counts_perchr.T
            counts_perchr[np.diag_indices_from(counts_perchr)] = counts_perchr[np.diag_indices_from(counts_perchr)] / 2
            counts_perchr = np.round(counts_perchr, 3)
            
            ## Output name for save
            if args.output is None:
                output_name = os.path.basename(args.filename)
                output_name = output_name.replace(".matrix", "_" + str(chrnames[i-1]) + "_dense.matrix")
            else:
                output_name = str(chrnames[i-1]) + "_" + args.output 
        
            if args.di is True:
                counts_perchr = np.hstack((bins[np.where(bins==str(chrnames[i-1]))[0]], counts_perchr))
            elif args.ins is True:
                ## add col header
                counts_perchr = np.hstack((bins_ins[np.where(bins==str(chrnames[i-1]))[0]], counts_perchr))
                ## add row header
                counts_perchr = np.vstack((np.insert(bins_ins[np.where(bins==str(chrnames[i-1]))[0]], 0, ''), counts_perchr))

            ## Save
            np.savetxt(output_name, counts_perchr, '%s', delimiter="\t")
