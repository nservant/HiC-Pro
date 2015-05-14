#!/usr/bin/python

# HiC-Pro
# Copyleft 2015 Institut Curie
# Author(s): Nicolas Servant, Eric Viara
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the
# GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to keep only valid 3C products - DE and SC are removed
"""

import getopt
import sys
import os
import re
import pysam
from bx.intervals.intersection import Intersecter, Interval


def usage():
    """Usage function"""
    print "Usage : python overlapMapped2HiCFragments.py"
    print "-f/--fragmentFile <Restriction fragment file GFF3>"
    print "-r/--mappedReadsFile <BAM/SAM file of mapped reads>"
    print "[-o/--outputDir] <Output directory. Default is current directory>"
    print "[-s/--shortestInsertSize] <Shortest insert size of mapped reads to consider>"
    print "[-l/--longestInsertSize] <Longest insert size of mapped reads to consider>"
    print "[-a/--all] <Write all additional output files, with information about the discarded reads>"
    print "[-S/--sam] <Output an additional SAM file with flag 'CT' for pairs classification>"
    print "[-v/--verbose] <Verbose>"
    print "[-h/--help] <Help>"
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "f:r:o:c:s:l:Svah",
            ["fragmentFile=",
             "mappedReadsFile=",
             "outputDir=", "cutSite=", "minInsertSize=",
             "maxInsertSize", "samOut", "verbose", "all", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts


def get_read_strand(read):
    """
    Conversion of read position to naive strand representation

    Parameters
    ----------
    read : list
        list of aligned reads
    """
    strand = "+"
    if read.is_reverse:
        strand = "-"
    return strand


def get_read_pos(read):
    """
    Return the read position (zero-based) used for the intersection with
    the restriction fragment

    The 5' end is not a good choice for the reverse reads (which contain part
    of the restriction site, and thus overlap the next restriction fragment)
    Using the left-most position (5' for forward, 3' for reverse) or the
    middle of the read should work but the middle of the reads might be more
    safe

    Parameters
    -----------
    read : list
        list of aligned reads
    """
    # 5 end
    # if (read.is_reverse):
    #    pos = read.pos + read.alen + 1 # (50 - 5) + 1 # zero-based transformation
    # else:
    #    pos = read.pos

    # Middle of the reads
    pos = read.pos + read.alen/2

    return pos


def get_ordered_reads(read1, read2):
    """
    Reorient reads

    The sequencing is usually not oriented. Reorient the reads so that r1 is
    always before r2

    read1 = [AlignedRead]
    read2 = [AlignedRead]
    """
    if read1.tid == read2.tid:
        if get_read_pos(read1) < get_read_pos(read2):
            r1 = read1
            r2 = read2
        else:
            r1 = read2
            r2 = read1
    else:
        if read1.tid < read2.tid:
            r1 = read1
            r2 = read2
        else:
            r1 = read2
            r2 = read1

    return r1, r2


def load_restriction_fragment(in_file, verbose):
    """
    Read a BED file and store the intervals in a tree

    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    in_file = input file [character]
    verbose = verbose mode [logical]

    """
    resFrag = {}
    if verbose:
        print "## Loading Restriction File Intervals '", in_file, "'..."

    bed_handle = open(in_file)
    for line in bed_handle:
        bedtab = line.split("\t")
        try:
            chromosome, start, end, name = bedtab[:4]
        except ValueError:
            # FIXME we might want a proper warning message here !
            continue

        # BED files are zero-based as Intervals objects
        start = int(start)  # + 1
        end = int(end)
        name = name.strip()
        if chromosome in resFrag.keys():
            tree = resFrag[chromosome]
            tree.add_interval(Interval(start, end, value={'name': name}))
        else:
            tree = Intersecter()
            tree.add_interval(Interval(start, end, value={'name': name}))
            resFrag[chromosome] = tree
    bed_handle.close()
    return resFrag


def get_overlapping_restriction_fragment(resFrag, chrom, read):
    """
    Intersect a given read with the set of restriction fragments

    ##
    resFrag = the restriction fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]

    """
    # Get 5' end
    pos = get_read_pos(read)
    
    if chrom in resFrag.keys():
        # Overlap with the 5' end of the read (zero-based)
        resfrag = resFrag[chrom].find(pos, pos+1)
        if len(resfrag) > 1:
            print "Error : ", len(resfrag),
            print " restriction fragments found for ", read.qname

        return resfrag[0]
    else:
        print "Warning - no restriction fragments for chromosome ", chrom
        return None


def is_self_circle(read1, read2):
    """
    Both reads are expected to be on the same restriction fragments

    Check the orientation of reads <-->
    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    ret = False
    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    # 1<- ->2 or 2<- ->1
    if get_read_strand(r1) == "-" and get_read_strand(r2) == "+":
        ret = True
    return ret


def is_dangling_end(read1, read2):
    """
    Both reads are expected to be on the same restriction fragments
    Check the orientation of reads -><-

    read1 : [AlignedRead]
    read2 : [AlignedRead]
    """
    ret = False
    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    # 1-> <-2 or 2-> <-1
    if get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        ret = True
    return ret


def get_valid_orientation(read1, read2):
    """
    Both reads are expected to be on the different restriction fragments

    Check the orientation of reads ->-> / <-<- / -><- / <-->

    read1 : [AlignedRead]
    read2 : [AlignedRead]

    """
    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)

    direction = None
    if get_read_strand(r1) == "+" and get_read_strand(r2) == "+":
        direction = "FF"
    elif get_read_strand(r1) == "-" and get_read_strand(r2) == "-":
        direction = "RR"
    elif get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        direction = "FR"
    elif get_read_strand(r1) == "-" and get_read_strand(r2) == "+":
        direction = "RF"

    return direction


def get_PE_fragment_size(read1, read2, resFrag1, resFrag2, interactionType):
    """
    Calculte the size of the DNA fragment library

    read1 : [AlignedRead]
    read2 : [AlignedRead]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    interactionType : Type of interaction from get_interaction_type() [str]

    """
    fragmentsize = None

    # Get oriented reads
    r1, r2 = get_ordered_reads(read1, read2)
    if not r1.is_unmapped and not r2.is_unmapped:
        if r1 == read2:
            rfrag1 = resFrag2
            rfrag2 = resFrag1
        else:
            rfrag1 = resFrag1
            rfrag2 = resFrag2

        r1pos = get_read_pos(r1)
        r2pos = get_read_pos(r2)

        if interactionType == "DE":
            fragmentsize = r2pos - r1pos
        elif interactionType == "SC":
            fragmentsize = (r1pos - rfrag1.start) + (rfrag2.end - r2pos)
        elif interactionType == "VI":
            if get_read_strand(r1) == "+":
                dr1 = rfrag1.end - r1pos
            else:
                dr1 = r1pos - rfrag1.start
            if get_read_strand(r2) == "+":
                dr2 = rfrag2.end - r2pos
            else:
                dr2 = r2pos - rfrag2.start
            fragmentsize = dr2 + dr1

    return fragmentsize


def get_interaction_type(read1, read1_chrom, resfrag1, read2,
                         read2_chrom, resfrag2, verbose):
    """
    Returns the interaction type

    For a given reads pair and their related restriction fragment, classify
    the 3C products as :

    - Interaction
    - Self circle
    - Dangling end
    - Unknown

    ##
    read1 = the R1 read of the pair [AlignedRead]
    read1_chrom = the chromosome of R1 read [character]
    resfrag1 = restrictin fragment overlapping the R1 read [interval]
    read2 = the R2 read of the pair [AlignedRead]
    read2_chrom = the chromosome of R2 read [character]
    resfrag2 = restrictin fragment overlapping the R2 read [interval]
    verbose = verbose mode [logical]

    """
    # If returned InteractionType=None -> Same restriction fragment
    # and same strand = Dump
    InteractionType = None

    if not (r1.is_unmapped) and not (r2.is_unmapped):
        # same restriction fragment
        if resfrag1 == resfrag2:
            # Self_circle <- ->
            if is_self_circle(read1, read2):
                InteractionType = "SC"
            # Dangling_end -> <-
            elif is_dangling_end(read1, read2):
                InteractionType = "DE"
        else:
            InteractionType = "VI"
    elif r1.is_unmapped or r2.is_unmapped:
        InteractionType = "SI"

    return InteractionType


def get_read_tag(read, tag):
    for t in read.tags:
        if t[0] == tag:
            return tag[1]
    return None


if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    inputFile = None
    outputFile = None
    samOut = False
    verbose = False
    allOutput = False
    minInsertSize = None
    maxInsertSize = None
    outputDir = "."

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--fragmentFile"):
            fragmentFile = arg
        elif opt in ("-r", "--mappedReadsFile"):
            mappedReadsFile = arg
        elif opt in ("-o", "--outputDir"):
            outputDir = arg
        elif opt in ("-s", "--shortestInsertSize"):
            minInsertSize = arg
        elif opt in ("-l", "--longestInsertSize"):
            maxInsertSize = arg
        elif opt in ("-a", "--all"):
            allOutput = True
        elif opt in ("-S", "--samOut"):
            samOut = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print "## overlapMapped2HiCFragments.py"
        print "## mappedReadsFile=", mappedReadsFile
        print "## fragmentFile=", fragmentFile
        print "## minInsertSize=", minInsertSize
        print "## maxInsertSize=", maxInsertSize
        print "## allOuput=", allOutput
        print "## SAM ouput=", samOut
        print "## verbose=", verbose, "\n"

    # Initialize variables
    reads_counter = 0
    de_counter = 0
    sc_counter = 0
    valid_counter = 0
    valid_counter_FF = 0
    valid_counter_RR = 0
    valid_counter_FR = 0
    valid_counter_RF = 0
    single_counter = 0
    dump_counter = 0

    baseReadsFile = os.path.basename(mappedReadsFile)
    baseReadsFile = re.sub(r'.bam|.sam', '', baseReadsFile)

    # Open handlers for output files
    handle_valid = open(outputDir + '/' + baseReadsFile + '.validPairs', 'w')

    if allOutput:
        handle_de = open(outputDir + '/' + baseReadsFile + '.DEPairs', 'w')
        handle_sc = open(outputDir + '/' + baseReadsFile + '.SCPairs', 'w')
        handle_dump = open(outputDir + '/' + baseReadsFile + '.DumpPairs', 'w')
        handle_single = open(outputDir + '/' + baseReadsFile + '.SinglePairs',
                             'w')

    # Read the BED file
    resFrag = load_restriction_fragment(fragmentFile, verbose)
    print resFrag.keys()
    
    # Read the SAM/BAM file
    if verbose:
        print "## Opening SAM/BAM file '", mappedReadsFile, "'..."
    samfile = pysam.Samfile(mappedReadsFile, "rb")

    if samOut:
        handle_sam = open(
            outputDir + '/' + baseReadsFile + '_interaction.sam', 'w')
        handle_sam = pysam.Samfile(
            outputDir + '/' + baseReadsFile + '_interaction.sam',
            "wh", header=samfile.header)

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    if verbose:
        print "## Classifying Interactions ..."

    for read in samfile.fetch(until_eof=True):
        reads_counter += 1
        cur_handler = None

        # First mate
        if read.is_read1:
            r1 = read
            if not r1.is_unmapped:
                r1_chrom = samfile.getrname(r1.tid)
                r1_resfrag = get_overlapping_restriction_fragment(resFrag, r1_chrom, r1)

        # Second mate
        elif read.is_read2:
            r2 = read
            if not r2.is_unmapped:
                r2_chrom = samfile.getrname(r2.tid)
                r2_resfrag = get_overlapping_restriction_fragment(resFrag, r2_chrom, r2)

            if r1_resfrag is not None and r2_resfrag is not None:
                interactionType = get_interaction_type(r1, r1_chrom, r1_resfrag, r2, r2_chrom, r2_resfrag, verbose)
                dist = get_PE_fragment_size(r1, r2, r1_resfrag, r2_resfrag, interactionType)

                # Check cut site in local mapping reads
                # if get_read_tag(r1, "RG") == "BML" or get_read_tag(r2, "RG") == "BML":
                #     bowloc_counter+=1
                #     if cutSite and (overlapRestrictionSite(r1, cutSite) or overlapRestrictionSite(r2, cutSite)):
                #         cutsite_counter+=1
                
                # Check Insert size criteria
                if (minInsertSize is not None and dist is not None and
                    dist < int(minInsertSize)) or \
                    (maxInsertSize is not None and dist is not None and dist > int(maxInsertSize)):
                    interactionType = "DUMP"

                if interactionType == "VI":
                    valid_counter += 1
                    cur_handler = handle_valid
                    validType = get_valid_orientation(r1, r2)
                    if validType == "RR":
                        valid_counter_RR += 1
                    elif validType == "FF":
                        valid_counter_FF += 1
                    elif validType == "FR":
                        valid_counter_FR += 1
                    elif validType == "RF":
                        valid_counter_RF += 1

                elif interactionType == "DE":
                    de_counter += 1
                    cur_handler = handle_de if allOutput else None

                elif interactionType == "SC":
                    sc_counter += 1
                    cur_handler = handle_sc if allOutput else None

                elif interactionType == "SI":
                    single_counter += 1
                    cur_handler = handle_single if allOutput else None
                else:
                    interactionType = "DUMP"
                    dump_counter += 1
                    cur_handler = handle_dump if allOutput else None
            else:
                interactionType = "DUMP"
                dump_counter += 1
                cur_handler = handle_dump if allOutput else None


            if cur_handler is not None:
                if not r1.is_unmapped and not r2.is_unmapped:
                    cur_handler.write(
                        r1.qname + "\t" +
                        r1_chrom + "\t" +
                        str(get_read_pos(r1)) + "\t" +
                        str(get_read_strand(r1)) + "\t" +
                        r2_chrom + "\t" +
                        str(get_read_pos(r2)) + "\t" +
                        str(get_read_strand(r2)) + "\t" +
                        str(dist) + "\n")
                elif r2.is_unmapped:
                    cur_handler.write(
                        r1.qname + "\t" +
                        r1_chrom + "\t" +
                        str(get_read_pos(r1)) + "\t" +
                        str(get_read_strand(r1)) + "\t" +
                        "+" + "\t" +
                        "0" + "\t" +
                        "*" + "\t" +
                        str(dist) + "\n")
                else:
                    cur_handler.write(
                        r1.qname + "\t" +
                        "*" + "\t" +
                        "0" + "\t" +
                        "*" + "\t" +
                        r2_chrom + "\t" +
                        str(get_read_pos(r2)) + "\t" +
                        str(get_read_strand(r2)) + "\t" +
                        str(dist) + "\n")

                if samOut:
                    r1.tags = r1.tags + [('CT', str(interactionType))]
                    r2.tags = r2.tags + [('CT', str(interactionType))]
                    handle_sam.write(r1)
                    handle_sam.write(r2)
                                   

            if (reads_counter % 100000 == 0 and verbose):
                print "##", reads_counter

    # Close handler
    handle_valid.close()
    if allOutput:
        handle_de.close()
        handle_sc.close()
        handle_dump.close()
        handle_single.close()

    # Write stats file
    handle_stat = open(outputDir + '/' + baseReadsFile + '.RSstat', 'w')
    handle_stat.write("## Hi-C processing\n")
    handle_stat.write("Valid_interaction_pairs\t" + str(valid_counter) + "\n")
    handle_stat.write(
        "Valid_interaction_pairs_FF\t" + str(valid_counter_FF) + "\n")
    handle_stat.write(
        "Valid_interaction_pairs_RR\t" + str(valid_counter_RR) + "\n")
    handle_stat.write(
        "Valid_interaction_pairs_RF\t" + str(valid_counter_RF) + "\n")
    handle_stat.write(
        "Valid_interaction_pairs_FR\t" + str(valid_counter_FR) + "\n")
    handle_stat.write("Dangling_end_pairs\t" + str(de_counter) + "\n")
    handle_stat.write("Self_Cycle_pairs\t" + str(sc_counter) + "\n")
    handle_stat.write("Single-end_pairs\t" + str(single_counter) + "\n")
    handle_stat.write("Dumped_pairs\t" + str(dump_counter) + "\n")
    handle_stat.close()

    if samOut:
        samfile.close()
