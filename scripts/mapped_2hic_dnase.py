#!/usr/bin/env python

# HiC-Pro
# Copyleft 2015 Institut Curie
# Author(s): Nicolas Servant, Eric Viara
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the
# GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to keep only valid pairs when no restriction enzyme are used (i.e. DNAse or Micro-HiC)
"""

import getopt
import sys
import os
import re
import pysam


def usage():
    """Usage function"""
    print("Usage : python mapped_2hic_dnase.py")
    print("-r/--mappedReadsFile <BAM/SAM file of mapped reads>")
    print("[-o/--outputDir] <Output directory. Default is current directory>")
    print("[-d/--minCisDist] <Minimum distance between intrachromosomal contact to consider>")
    print("[-g/--gtag] <Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification>")
    print("[-a/--all] <Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "r:o:d:g:avh",
            ["mappedReadsFile=",
             "outputDir=", "minDist=", "gatg", "all", "verbose", "help"])
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


def get_read_pos(read, st="start"):
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
    if st == "middle":
        pos = read.reference_start + int(read.alen/2)
    elif st =="start":
        pos = get_read_start(read)
    elif st == "left":
        pos = read.reference_start

    return pos


def get_read_start(read):
    """                                                                                                                                                                                                        
    Return the 5' end of the read                                                                                                                                                                              
    """
    if read.is_reverse:
        pos = read.reference_start + read.alen -1
    else:
        pos = read.reference_start
    return pos


def get_ordered_reads(read1, read2):
    """
    Reorient reads

    The sequencing is usually not oriented. Reorient the reads so that r1 is
    always before r2

    read1 = [AlignedRead]
    read2 = [AlignedRead]
    """
    if read1.reference_id == read2.reference_id:
        if get_read_pos(read1) < get_read_pos(read2):
            r1, r2 = read1, read2
        else:
            r1, r2 = read2, read1
    else:
        if read1.reference_id < read2.reference_id:
            r1, r2 = read1, read2
        else:
            r1, r2 = read2, read1

    return r1, r2


def isIntraChrom(read1, read2):
    """
    Return true is the reads pair is intrachromosomal
    
    read1 : [AlignedRead]
    read2 : [AlignedRead]

    """
    if read1.reference_id == read2.reference_id:
        return True
    else:
        return False


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


def get_cis_dist(read1, read2):
     """
     Calculte the size of the DNA fragment library

     read1 : [AlignedRead]
     read2 : [AlignedRead]

     """
     # Get oriented reads
     ##r1, r2 = get_ordered_reads(read1, read2)
     dist = None
     if not r1.is_unmapped and not r2.is_unmapped:         
         ## Contact distances can be calculated for intrachromosomal reads only
         if isIntraChrom(read1, read2):
             r1pos = get_read_pos(read1)
             r2pos = get_read_pos(read2)
             dist = abs(r1pos - r2pos)
     return dist


def get_read_tag(read, tag):
    for t in read.get_tags():
        if t[0] == tag:
            return t[1]
    return None


if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    verbose = False
    allOutput = False
    minInsertSize = None
    maxInsertSize = None
    minDist = None
    outputDir = "."
    gtag = None

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-r", "--mappedReadsFile"):
            mappedReadsFile = arg
        elif opt in ("-o", "--outputDir"):
            outputDir = arg
        elif opt in ("-d", "--minCisDist"):
            minDist = arg
        elif opt in ("-g", "--gtag"):
            gtag = arg
        elif opt in ("-a", "--all"):
            allOutput = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print("## overlapMapped2HiCFragments.py")
        print("## mappedReadsFile=", mappedReadsFile)
        print("## minCisDist=", minDist)
        print("## allOuput=", allOutput)
        print("## verbose={}\n".format(verbose))

    # Initialize variables
    reads_counter = 0
    valid_counter = 0
    valid_counter_FF = 0
    valid_counter_RR = 0
    valid_counter_FR = 0
    valid_counter_RF = 0
    single_counter = 0
    dump_counter = 0
    filt_counter = 0

    # AS counter
    G1G1_ascounter = 0
    G2G2_ascounter = 0
    G1U_ascounter = 0
    UG1_ascounter = 0
    G2U_ascounter = 0
    UG2_ascounter = 0
    G1G2_ascounter = 0
    G2G1_ascounter = 0
    UU_ascounter = 0
    CF_ascounter = 0

    baseReadsFile = os.path.basename(mappedReadsFile)
    baseReadsFile = re.sub(r'\.bam$|\.sam$', '', baseReadsFile)

    # Open handlers for output files
    handle_valid = open(outputDir + '/' + baseReadsFile + '.validPairs', 'w')

    if allOutput:
        handle_dump = open(outputDir + '/' + baseReadsFile + '.DumpPairs', 'w')
        handle_single = open(outputDir + '/' + baseReadsFile + '.SinglePairs','w')
        handle_filt = open(outputDir + '/' + baseReadsFile + '.FiltPairs','w')

    # Read the SAM/BAM file
    if verbose:
        print("## Opening SAM/BAM file {} ...".format(mappedReadsFile))
    samfile = pysam.Samfile(mappedReadsFile, "rb")

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    for read in samfile.fetch(until_eof=True):
        reads_counter += 1
        cur_handler = None
        interactionType = None
        htag = ""

        # First mate
        if read.is_read1:
            r1 = read
            if not r1.is_unmapped:
                r1_chrom = samfile.get_reference_name(r1.reference_id)
            else:
                r1_chrom = None

        # Second mate
        elif read.is_read2:
            r2 = read
            if not r2.is_unmapped:
                r2_chrom = samfile.get_reference_name(r2.reference_id)
            else:
                r2_chrom = None

            if isIntraChrom(r1, r2):
                dist = get_cis_dist(r1, r2)
            else:
                dist = None

            # Check singleton
            if r1.is_unmapped or r2.is_unmapped:
                interactionType = "SI"
                single_counter += 1
                cur_handler = handle_single if allOutput else None

            # Check Distance criteria - Filter
            if (minDist is not None and dist is not None and dist < int(minDist)):
                interactionType = "FILT"
                filt_counter += 1
                cur_handler = handle_filt if allOutput else None

            # By default pair is valid
            if interactionType == None:
                interactionType = "VI"
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
                else:
                    interactionType = "DUMP"
                    dump_counter += 1
                    cur_handler = handle_dump if allOutput else None



            # Split valid pairs based on XA tag
            if gtag is not None:
                r1as = get_read_tag(r1, gtag)
                r2as = get_read_tag(r2, gtag)
                        
                if r1as == 1 and r2as == 1:
                    G1G1_ascounter += 1
                elif r1as == 2 and r2as == 2:
                    G2G2_ascounter += 1
                elif r1as == 1 and r2as == 0:
                    G1U_ascounter += 1
                elif r1as == 0 and r2as == 1:
                    UG1_ascounter += 1
                elif r1as == 2 and r2as == 0:
                    G2U_ascounter += 1
                elif r1as == 0 and r2as == 2:
                    UG2_ascounter += 1
                elif r1as == 1 and r2as == 2:
                    G1G2_ascounter += 1
                elif r1as == 2 and r2as == 1:
                    G2G1_ascounter += 1
                elif r1as == 3 or r2as == 3:
                    CF_ascounter += 1
                else:
                    UU_ascounter += 1
                        
       
            if cur_handler is not None:
                if not r1.is_unmapped and not r2.is_unmapped:
                    
                    ##reorient reads to ease duplicates removal
                    or1, or2 = get_ordered_reads(r1, r2)
                    or1_chrom = samfile.get_reference_name(or1.reference_id)
                    or2_chrom = samfile.get_reference_name(or2.reference_id)

                    ##reset as tag now that the reads are oriented
                    r1as = get_read_tag(or1, gtag)
                    r2as = get_read_tag(or2, gtag)
                    if gtag is not None:
                        htag = str(r1as)+"-"+str(r2as)
                        
                    cur_handler.write(
                        or1.query_name + "\t" +
                        or1_chrom + "\t" +
                        str(get_read_pos(or1)+1) + "\t" +
                        str(get_read_strand(or1)) + "\t" +
                        or2_chrom + "\t" +
                        str(get_read_pos(or2)+1) + "\t" +
                        str(get_read_strand(or2)) + "\t" +
                        "NA" + "\t" + ##dist 
                        "NA" + "\t" + ##resfrag1
                        "NA" + "\t" + ##resfrag2
                        str(or1.mapping_quality) + "\t" + 
                        str(or2.mapping_quality) + "\t" + 
                        str(htag) + "\n")
                
                elif r2.is_unmapped and not r1.is_unmapped:
                    cur_handler.write(
                        r1.query_name + "\t" +
                        r1_chrom + "\t" +
                        str(get_read_pos(r1)+1) + "\t" +
                        str(get_read_strand(r1)) + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" + 
                        "*" + "\t" +
                        "*" + "\t" +
                        str(r1.mapping_quality) + "\t" + 
                        "*" + "\n")
                elif r1.is_unmapped and not r2.is_unmapped:
                    cur_handler.write(
                        r2.query_name + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        r2_chrom + "\t" +
                        str(get_read_pos(r2)+1) + "\t" +
                        str(get_read_strand(r2)) + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" + 
                        str(r2.mapping_quality) + "\n")

            if (reads_counter % 100000 == 0 and verbose):
                print("##", reads_counter)

    # Close handler
    handle_valid.close()
    if allOutput:
        handle_dump.close()
        handle_single.close()
        handle_filt.close()

    # Write stats file
    with open(outputDir + '/' + baseReadsFile + '.RSstat', 'w') as handle_stat:
        handle_stat.write("## Hi-C processing - no restriction fragments\n")
        handle_stat.write("Valid_interaction_pairs\t" + str(valid_counter) + "\n")
        handle_stat.write("Valid_interaction_pairs_FF\t" + str(valid_counter_FF) + "\n")
        handle_stat.write("Valid_interaction_pairs_RR\t" + str(valid_counter_RR) + "\n")
        handle_stat.write("Valid_interaction_pairs_RF\t" + str(valid_counter_RF) + "\n")
        handle_stat.write("Valid_interaction_pairs_FR\t" + str(valid_counter_FR) + "\n")
        handle_stat.write("Single-end_pairs\t" + str(single_counter) + "\n")
        handle_stat.write("Filtered_pairs\t" + str(filt_counter) + "\n")
        handle_stat.write("Dumped_pairs\t" + str(dump_counter) + "\n")

    ## Write AS report
        if gtag is not None:
            handle_stat.write("## ======================================\n")
            handle_stat.write("## Allele specific information\n")
            handle_stat.write("Valid_pairs_from_ref_genome_(1-1)\t" + str(G1G1_ascounter) + "\n")
            handle_stat.write("Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)\t" + str(UG1_ascounter+G1U_ascounter) + "\n")
            handle_stat.write("Valid_pairs_from_alt_genome_(2-2)\t" + str(G2G2_ascounter) + "\n")
            handle_stat.write("Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)\t" + str(UG2_ascounter+G2U_ascounter) + "\n")
            handle_stat.write("Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)\t" + str(G1G2_ascounter+G2G1_ascounter) + "\n")
            handle_stat.write("Valid_pairs_with_both_unassigned_mated_(0-0)\t" + str(UU_ascounter) + "\n")
            handle_stat.write("Valid_pairs_with_at_least_one_conflicting_mate_(3-)\t" + str(CF_ascounter) + "\n")



