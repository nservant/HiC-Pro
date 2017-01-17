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
    infile = open(filename)
    chr_id = None
    big_str = ""
    indices = []
    all_indices = []
    contig_names = []
    c = 0
    for line in infile:
        c += 1
        if line.startswith(">"):
            print line.split()[0][1:], "..."
            # If this is not the first chromosome, find the indices and append
            # them to the list
            if chr_id is not None:
                for rs in range(len(sequences)):
                    pattern = "(?=%s)" % sequences[rs].lower()
                    indices += [m.start() + offset[rs]
                                for m in re.finditer(pattern, big_str)]
                indices.sort()
                all_indices.append(indices)
                indices = []
            # This is a new chromosome. Empty the sequence string, and add the
            # correct chrom id
            big_str = ""
            chr_id = line.split()[0][1:]
            if chr_id in contig_names:
                print "The fasta file contains several instance of",
                print chr_id, ". Exit."
                sys.exit(-1)
            contig_names.append(chr_id)
        else:
            # As long as we don't change chromosomes, continue reading the
            # file, and appending the sequences
            big_str += line.lower().strip()
    # Add the indices for the last chromosome
    for rs in range(len(sequences)):
        pattern = "(?=%s)" % sequences[rs].lower()
        indices += [m.start() + offset[rs]
                    for m in re.finditer(pattern, big_str)]
    indices.sort()
    all_indices.append(indices)
    return contig_names, all_indices


def find_chromsomose_lengths(reference_filename):
    chromosome_lengths = []
    chromosome_names = []
    length = None
    infile = open(reference_filename)
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
    cutsites = args.res_sites

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
            print "Unable to detect offset for", cseq
            print "Please, use '^' to specified the cutting position,",
            print "i.e A^GATCT for HindIII digestion"
            sys.exit(-1)
        offset.append(offpos)
        sequences.append(re.sub('\^', '', cseq))

    if out is None:
        out = os.path.splitext(filename)[0] + "_fragments.bed"

    print "Analyzing", filename
    print "Restriction site(s)", ",".join(sequences)
    print "Offset(s)",  ','.join(str(x) for x in offset)

    # Read fasta file and look for rs per chromosome
    contig_names, all_indices = find_re_sites(filename, sequences,
                                              offset=offset)
    _, lengths = find_chromsomose_lengths(filename)

    valid_fragments = []
    for i, indices in enumerate(all_indices):
        valid_fragments_chr = np.concatenate(
            [np.concatenate([[0], indices])[:, np.newaxis],
             np.concatenate([indices, [lengths[i]]])[:, np.newaxis]],
            axis=1)
        valid_fragments.append(valid_fragments_chr)

    # Write results
    print "Writing to", out, "..."
    outfile = open(out, "w")
    for chrom_name, indices in zip(contig_names, valid_fragments):
        frag_id = 0
        for begin, end in indices:
            # allow to remove cases where the enzyme cut at
            # the first position of the chromosome
            if end > begin:
                frag_id += 1
                frag_name = "HIC_%s_%d" % (chrom_name, frag_id)
                outfile.write(
                    "%s\t%d\t%d\t%s\t0\t+\n" % (chrom_name, begin,
                                                end, frag_name))
    outfile.close()
