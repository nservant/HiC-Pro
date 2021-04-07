#!/usr/bin/env python

# HiC-Pro
# Copyleft 2015 Institut Curie
# Author(s): Nelle Varoquaux, Nicolas Servant
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the
# GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to extract restriction fragment from a fasta file and output a BED file
"""

import argparse
import re
import os
import sys
import numpy as np

RE_cutsite = {
    "mboi": ["^GATC"],
    "dpnii": ["^GATC"],
    "bglii": ["A^GATCT"],
    "hindiii": ["A^AGCTT"]}


def find_re_sites(filename, sequences, offset):
    with open(filename, 'r') as infile:
        chr_id = None
        big_str = ""
        indices = []
        all_indices = []
        contig_names = []
        c = 0
        for line in infile:
            c += 1
            if line.startswith(">"):
                print("{}...".format(line.split()[0][1:]))
                # If this is not the first chromosome, find the indices and append
                # them to the list
                if chr_id is not None:
                     for rs in range(len(sequences)):
                         pattern = "(?={})".format(sequences[rs].lower())
                         indices += [m.start() + offset[rs]\
                         for m in re.finditer(pattern, big_str)]
                     indices.sort()
                     all_indices.append(indices)
                     indices = []

                # This is a new chromosome. Empty the sequence string, and add the
                # correct chrom id
                big_str = ""
                chr_id = line.split()[0][1:]
                if chr_id in contig_names:
                    print("The fasta file contains several instance of {}. Exit.".format(chr_id))
                    sys.exit(-1)
                contig_names.append(chr_id)
            else:
                # As long as we don't change chromosomes, continue reading the
                # file, and appending the sequences
                big_str += line.lower().strip()
        # Add the indices for the last chromosome
        for rs in range(len(sequences)):
            pattern = "(?={})".format(sequences[rs].lower())
            indices += [m.start() + offset[rs]
                        for m in re.finditer(pattern, big_str)]
        indices.sort()
        all_indices.append(indices)
    
    return contig_names, all_indices


def find_chromsomose_lengths(reference_filename):
    chromosome_lengths = []
    chromosome_names = []
    length = None
    with open(reference_filename, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                chromosome_names.append(line[1:].strip())
                if length is not None:
                    chromosome_lengths.append(length)
                length = 0
            else:
                length += len(line.strip())
        chromosome_lengths.append(length)
    return chromosome_names, np.array(chromosome_lengths)


def replaceN(cs):
    npos = int(cs.find('N'))
    cseql = []
    if npos != -1:
        for nuc in ["A","C","G","T"]:
            tmp = cs.replace('N', nuc, 1)
            tmpl = replaceN(tmp)
            if type(tmpl) == list:
                cseql = cseql + tmpl
            else:
                cseql.append(tmpl)
    else:
        cseql.append(cs)
    return cseql


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile')
    parser.add_argument('-r', '--restriction_sites',
                        dest='res_sites',
                        nargs='+',
                        help=("The cutting position has to be specified using "
                              "'^'. For instance, -r A^AGCTT for HindIII "
                              "digestion. Several restriction enzyme can be "
                              "specified."))
    parser.add_argument('-o', '--out', default=None)
    args = parser.parse_args()

    filename = args.fastafile
    out = args.out
    
    # Split restriction sites if comma-separated
    cutsites=[]
    for s in args.res_sites:
        for m in s.split(','):
            cutsites.append(m)
                
    # process args and get restriction enzyme sequences
    sequences = []
    offset = []
    for cs in cutsites:
        if cs.lower() in RE_cutsite:
            cseq = ''.join(RE_cutsite[cs.lower()])
        else:
            cseq = cs

        offpos = int(cseq.find('^'))
        if offpos == -1:
            print("Unable to detect offset for {}. Please, use '^' to specify the cutting position,\
                   i.e A^GATCT for HindIII digestion.".format(cseq))
            sys.exit(-1)

        for nuc in list(set(cseq)):
            if nuc not in ['A','T','G','C','N','^']:
                print("Find unexpected character ['{}']in restriction motif".format(nuc))
                print("Note that multiple motifs should be separated by a space (not a comma !)")

                sys.exit(-1)

        offset.append(offpos)
        sequences.append(re.sub('\^', '', cseq))

    # replace all N in restriction motif
    sequences_without_N = []
    offset_without_N = []
    for rs in range(len(sequences)):
        nrs = replaceN(sequences[rs])
        sequences_without_N = sequences_without_N + nrs
        offset_without_N = offset_without_N + [offset[rs]] * len(nrs)
          
    sequences = sequences_without_N
    offset = offset_without_N
    
    if out is None:
        out = os.path.splitext(filename)[0] + "_fragments.bed"

    print("Analyzing", filename)
    print("Restriction site(s)", ",".join(sequences))
    print("Offset(s)",  ','.join(str(x) for x in offset))

    # Read fasta file and look for rs per chromosome
    contig_names, all_indices = find_re_sites(filename, sequences,  offset=offset)
    _, lengths = find_chromsomose_lengths(filename)

    valid_fragments = []
    for i, indices in enumerate(all_indices):
        valid_fragments_chr = np.concatenate(
            [np.concatenate([[0], indices])[:, np.newaxis],
             np.concatenate([indices, [lengths[i]]])[:, np.newaxis]],
            axis=1)
        valid_fragments.append(valid_fragments_chr)

    # Write results
    print("Writing to {} ...".format(out))
    with open(out, 'w') as outfile:
        for chrom_name, indices in zip(contig_names, valid_fragments):
            frag_id = 0
            for begin, end in indices:
                # allow to remove cases where the enzyme cut at
                # the first position of the chromosome
                if end > begin:
                    frag_id += 1
                    frag_name = "HIC_{}_{}".format(str(chrom_name), int(frag_id))
                    outfile.write("{}\t{}\t{}\t{}\t0\t+\n".format(str(chrom_name), int(begin), int(end), str(frag_name)))
