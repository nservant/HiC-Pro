#!/usr/bin/env python

# HiC-Pro
# Copyleft 2015 Institut Curie
# Author(s): Nicolas Servant, Eric Viara
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the
# GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to keep only valid 3C products - DE and SC are removed
Output is : readname / 
"""
import time
import getopt
import sys
import os
import re
import pysam
from bx.intervals.intersection import Intersecter, Interval


def usage():
    """Usage function"""
    print("Usage : python mapped_2hic_fragments.py")
    print("-f/--fragmentFile <Restriction fragment file GFF3>")
    print("-r/--mappedReadsFile <BAM/SAM file of mapped reads>")
    print("[-o/--outputDir] <Output directory. Default is current directory>")
    print("[-s/--shortestInsertSize] <Shortest insert size of mapped reads to consider>")
    print("[-l/--longestInsertSize] <Longest insert size of mapped reads to consider>")
    print("[-t/--shortestFragmentLength] <Shortest restriction fragment length to consider>")
    print("[-m/--longestFragmentLength] <Longest restriction fragment length to consider>")
    print("[-d/--minCisDist] <Minimum distance between intrachromosomal contact to consider>")
    print("[-g/--gtag] <Genotype tag. If specified, this tag will be reported in the valid pairs output for allele specific classification>")
    print("[-a/--all] <Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)>")
    print("[-S/--sam] <Output an additional SAM file with flag 'CT' for pairs classification>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "f:r:o:s:l:t:m:d:g:Svah",
            ["fragmentFile=",
             "mappedReadsFile=",
             "outputDir=", 
             "minInsertSize=", "maxInsertSize", 
             "minFragSize", "maxFragSize", 
             "minDist",
             "gatg", "sam", "verbose", "all", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts


def timing(function, *args):
    """
    Run a fonction and eturn the run time and the result of the function
    If the function requires arguments, those can be passed in
    """
    startTime = time.time()
    result = function(*args)
    print('{} function took {:.3f}ms'.format(function.__name__, (time.time() - startTime) * 1000))
    return result


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


def isIntraChrom(read1, read2):
    """
    Return true is the reads pair is intrachromosomal
    
    read1 : [AlignedRead]
    read2 : [AlignedRead]

    """
    if read1.tid == read2.tid:
        return True
    return False


def get_cis_dist(read1, read2):
     """
     Calculte the contact distance between two intrachromosomal reads

     read1 : [AlignedRead]
     read2 : [AlignedRead]

     """
     # Get oriented reads
     ##r1, r2 = get_ordered_reads(read1, read2)
     dist = None
     if not read1.is_unmapped and not read2.is_unmapped:         
         ## Contact distances can be calculated for intrachromosomal reads only
         if isIntraChrom(read1, read2):
             r1pos, r2pos = get_read_pos(read1), get_read_pos(read2)
             dist = abs(r1pos - r2pos)
     return dist


def get_read_pos(read, st="start"):
    """
    Return the read position (zero-based) used for the intersection with
    the restriction fragment

    The 5' end is not a good choice for the reverse reads (which contain part
    of the restriction site, and thus overlap the next restriction fragment)
    Using the left-most position (ie. start, 5' for forward, 3' for reverse) or the
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
    always before r2.
    Sequencing is always performed from 5' to 3' end
    So in unstranded case, we can have

    1              2
    --->           --->
    ==========  or =========
         <----          <---
             2             1

    Reordering the reads allow to always be in the first case
    read1 = [AlignedRead]
    read2 = [AlignedRead]
    """
    if read1.tid == read2.tid:
        if get_read_pos(read1) < get_read_pos(read2):
            r1, r2 = read1, read2
        else:
            r1, r2 = read2, read1
    else:
        if read1.tid < read2.tid:
            r1, r2 = read1, read2
        else:
            r1, r2 = read2, read1
                
    return r1, r2

def load_restriction_fragment(in_file, minfragsize=None, maxfragsize=None, verbose=False):
    """
    Read a BED file and store the intervals in a tree

    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    in_file = input file [character]
    verbose = verbose mode [logical]

    """
    resFrag = {}
    if verbose:
        print("## Loading Restriction File Intervals {} ...".format(in_file))
    bed_handle = open(in_file)
    nline = 0
    nfilt = 0
    for line in bed_handle:
         nline += 1
         bedtab = line.split("\t")
         try:
              chromosome, start, end, name = bedtab[:4]
         except ValueError:
              print("Warning : wrong input format in line {}. Not a BED file ?!".format(nline))
              continue

        # BED files are zero-based as Intervals objects
         start = int(start)  # + 1
         end = int(end)
         fragl = abs(end - start)
         name = name.strip()

         ## Discard fragments outside the size range
         filt = False
         if minfragsize != None and int(fragl) < int(minfragsize):
             nfilt += 1
             filt = True
         elif maxfragsize != None and int(fragl) > int(maxfragsize):
             nfilt += 1
             filt = True
       
         if chromosome in resFrag:
             tree = resFrag[chromosome]
             tree.add_interval(Interval(start, end, value={'name': name, 'filter': filt}))
         else:
             tree = Intersecter()
             tree.add_interval(Interval(start, end, value={'name': name, 'filter': filt}))
             resFrag[chromosome] = tree
    
    if nfilt > 0:
        print("Warning : {} fragment(s) outside of range and discarded. {} remaining.".format(nfilt, nline - nfilt))
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
    # Get read position (middle or start)
    pos = get_read_pos(read, st="middle")
    
    if chrom in resFrag:
        # Overlap with the position of the read (zero-based)
        resfrag = resFrag[chrom].find(pos, pos+1)
        if len(resfrag) > 1:
            print("Warning : {} restictions fragments found for {} -skipped".format(len(resfrag), read.query_name))
            return None
        elif len(resfrag) == 0:
            print("Warning - no restriction fragments for {} at {} : {}".format(read.query_name, chrom, pos))
            return None
        else:
            return resfrag[0]
    else:
        print("Warning - no restriction fragments for {} at {} : {}".format(read.qname, chrom, pos))
        return None


def are_contiguous_fragments(frag1, frag2, chr1, chr2):
    '''
    Compare fragment positions to check if they are contiguous
    '''
    ret = False
    if chr1 == chr2:
        if int(frag1.start) < int(frag2.start):
            d = int(frag2.start) - int(frag1.end)
        else:
            d = int(frag1.start) - int(frag2.end)
            
        if d == 0:
            ret = True
    
    return ret

def is_religation(read1, read2, frag1, frag2):
    """
    Reads are expected to map adjacent fragments
    Check the orientation of reads -><-

    """
    ret = False
    if are_contiguous_fragments(frag1, frag2, read1.tid, read2.tid):
        #r1, r2 = get_ordered_reads(read1, read2)
        #if get_read_strand(r1) == "+" and get_read_strand(r2) == "-":
        ret = True
    return ret


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
    resfrag1 = restriction fragment overlapping the R1 read [interval]
    resfrag1 = restriction fragment overlapping the R1 read [interval]
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

        ## In this case use the read start !
        r1pos = get_read_start(r1)
        r2pos = get_read_start(r2)

        if interactionType == "DE" or interactionType == "RE":
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
    - Religation
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
    interactionType = None
      
    if not read1.is_unmapped and not read2.is_unmapped and resfrag1 is not None and resfrag2 is not None:
        # same restriction fragment
        if resfrag1 == resfrag2:
            # Self_circle <- ->
            if is_self_circle(read1, read2):
                interactionType = "SC"
            # Dangling_end -> <-
            elif is_dangling_end(read1, read2):
                interactionType = "DE"
        elif is_religation(read1, read2, resfrag1, resfrag2):
            interactionType = "RE"
        else:
            interactionType = "VI"
    elif r1.is_unmapped or r2.is_unmapped:
        interactionType = "SI"

    return interactionType


def get_read_tag(read, tag):
    for t in read.get_tags():
        if t[0] == tag:
            return t[1]
    return None


if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    samOut = False
    verbose = False
    allOutput = False
    minInsertSize = None
    maxInsertSize = None
    minFragSize = None
    maxFragSize = None
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
        elif opt in ("-t", "--shortestFragmentLength"):
            minFragSize = arg
        elif opt in ("-m", "--longestFragmentLength"):
            maxFragSize = arg
        elif opt in ("-d", "--minCisDist"):
            minDist = arg
        elif opt in ("-g", "--gtag"):
            gtag = arg
        elif opt in ("-a", "--all"):
            allOutput = True
        elif opt in ("-S", "--sam"):
            samOut = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print("## overlapMapped2HiCFragments.py")
        print("## mappedReadsFile=", mappedReadsFile)
        print("## fragmentFile=", fragmentFile)
        print("## minInsertSize=", minInsertSize)
        print("## maxInsertSize=", maxInsertSize)
        print("## minFragSize=", minFragSize)
        print("## maxFragSize=", maxFragSize)
        print("## allOuput=", allOutput)
        print("## SAM ouput=", samOut)
        print("## verbose={}\n".format(verbose))

    # Initialize variables
    reads_counter = 0
    de_counter = 0
    re_counter = 0
    sc_counter = 0
    valid_counter = 0
    valid_counter_FF = 0
    valid_counter_RR = 0
    valid_counter_FR = 0
    valid_counter_RF = 0
    single_counter = 0
    dump_counter = 0
    filt_counter = 0

    ## AS counter
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
        handle_de = open(outputDir + '/' + baseReadsFile + '.DEPairs', 'w')
        handle_re = open(outputDir + '/' + baseReadsFile + '.REPairs', 'w')
        handle_sc = open(outputDir + '/' + baseReadsFile + '.SCPairs', 'w')
        handle_dump = open(outputDir + '/' + baseReadsFile + '.DumpPairs', 'w')
        handle_single = open(outputDir + '/' + baseReadsFile + '.SinglePairs', 'w')
        handle_filt = open(outputDir + '/' + baseReadsFile + '.FiltPairs', 'w')

    # Read the BED file
    resFrag = timing(load_restriction_fragment, fragmentFile, minFragSize, maxFragSize, verbose)
     
    # Read the SAM/BAM file
    if verbose:
        print("## Opening SAM/BAM file {} ...".format(mappedReadsFile))
    samfile = pysam.Samfile(mappedReadsFile, "rb")

    if samOut:
        handle_sam = pysam.AlignmentFile(outputDir + '/' + baseReadsFile + '_interaction.bam', "wb", template=samfile)

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    if verbose:
        print("## Classifying Interactions ...")

    for read in samfile.fetch(until_eof=True):
        reads_counter += 1
        cur_handler = None
        htag = ""

        # First mate
        if read.is_read1:
            r1 = read
            if not r1.is_unmapped:
                r1_chrom = samfile.get_reference_name(r1.tid)
                r1_resfrag = get_overlapping_restriction_fragment(resFrag, r1_chrom, r1)
            else:
                r1_resfrag = None
                r1_chrom = None

        # Second mate
        elif read.is_read2:
            r2 = read
            if not r2.is_unmapped:
                r2_chrom = samfile.get_reference_name(r2.tid)
                r2_resfrag = get_overlapping_restriction_fragment(resFrag, r2_chrom, r2)
            else:
                r2_resfrag = None
                r2_chrom = None

            if r1_resfrag is not None or r2_resfrag is not None:

                interactionType = get_interaction_type(r1, r1_chrom, r1_resfrag, r2, r2_chrom, r2_resfrag, verbose)
                dist = get_PE_fragment_size(r1, r2, r1_resfrag, r2_resfrag, interactionType)
                cdist = get_cis_dist(r1, r2)
                
                ## Filter based on restriction fragments
                if (r1_resfrag is not None and r1_resfrag.value['filter'] == True) or (r2_resfrag is not None and r2_resfrag.value['filter']) == True:
                    interactionType = "FILT"
   
                # Check Insert size criteria - FILT
                if (minInsertSize is not None and dist is not None and
                    dist < int(minInsertSize)) or \
                    (maxInsertSize is not None and dist is not None and dist > int(maxInsertSize)):
                    interactionType = "FILT"

                # Check Distance criteria - FILT
                # Done for VI otherwise this criteria will overwrite all other invalid classification
                if (interactionType == "VI" and minDist is not None and cdist is not None and cdist < int(minDist)):
                    interactionType = "FILT"
        
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

                    ## Counts valid pairs based on XA tag
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

                elif interactionType == "DE":
                    de_counter += 1
                    cur_handler = handle_de if allOutput else None

                elif interactionType == "RE":
                    re_counter += 1
                    cur_handler = handle_re if allOutput else None

                elif interactionType == "SC":
                    sc_counter += 1
                    cur_handler = handle_sc if allOutput else None

                elif interactionType == "SI":
                    single_counter += 1
                    cur_handler = handle_single if allOutput else None
                
                elif interactionType == "FILT":
                    filt_counter += 1
                    cur_handler = handle_filt if allOutput else None
                
                else:
                    interactionType = "DUMP"
                    dump_counter += 1
                    cur_handler = handle_dump if allOutput else None
            else:
                interactionType = "DUMP"
                dump_counter += 1
                cur_handler = handle_dump if allOutput else None
                dist = None

            ## Write results in right handler
            if cur_handler is not None:
                if not r1.is_unmapped and not r2.is_unmapped:                 
                    ##reorient reads to ease duplicates removal
                    or1, or2 = get_ordered_reads(r1, r2)
                    or1_chrom = samfile.get_reference_name(or1.tid)
                    or2_chrom = samfile.get_reference_name(or2.tid)
                    
                    ##reset as tag now that the reads are oriented
                    r1as = get_read_tag(or1, gtag)
                    r2as = get_read_tag(or2, gtag)
                    if gtag is not None:
                        htag = str(r1as)+"-"+str(r2as)

                    ##get fragment name and reorient if necessary
                    if or1 == r1 and or2 == r2:
                        or1_resfrag = r1_resfrag
                        or2_resfrag = r2_resfrag
                    elif or1 == r2 and or2 == r1:
                        or1_resfrag = r2_resfrag
                        or2_resfrag = r1_resfrag

                    if or1_resfrag is not None:
                        or1_fragname = or1_resfrag.value['name']
                    else:
                        or1_fragname = 'None'
                        
                    if or2_resfrag is not None:
                        or2_fragname = or2_resfrag.value['name']
                    else:
                        or2_fragname = 'None'
                        
                    cur_handler.write(
                        or1.query_name + "\t" +
                        or1_chrom + "\t" +
                        str(get_read_pos(or1)+1) + "\t" +
                        str(get_read_strand(or1)) + "\t" +
                        or2_chrom + "\t" +
                        str(get_read_pos(or2)+1) + "\t" +
                        str(get_read_strand(or2)) + "\t" +
                        str(dist) + "\t" + 
                        or1_fragname + "\t" +
                        or2_fragname + "\t" +
                        str(or1.mapping_quality) + "\t" + 
                        str(or2.mapping_quality) + "\t" + 
                        str(htag) + "\n")

                elif r2.is_unmapped and not r1.is_unmapped:
                    if r1_resfrag is not None:
                        r1_fragname = r1_resfrag.value['name']
                          
                    cur_handler.write(
                        r1.query_name + "\t" +
                        r1_chrom + "\t" +
                        str(get_read_pos(r1)+1) + "\t" +
                        str(get_read_strand(r1)) + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" +
                        "*" + "\t" + 
                        r1_fragname + "\t" +
                        "*" + "\t" +
                        str(r1.mapping_quality) + "\t" + 
                        "*" + "\n")
                elif r1.is_unmapped and not r2.is_unmapped:
                    if r2_resfrag is not None:
                        r2_fragname = r2_resfrag.value['name']
                    
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
                        r2_fragname + "\t" +
                        "*" + "\t" +
                        str(r2.mapping_quality) + "\n")

                ## Keep initial order    
                if samOut:
                    r1.tags = r1.tags + [('CT', str(interactionType))]
                    r2.tags = r2.tags + [('CT', str(interactionType))]
                    handle_sam.write(r1)
                    handle_sam.write(r2)

            if (reads_counter % 100000 == 0 and verbose):
                print("##", reads_counter)

    # Close handler
    handle_valid.close()
    if allOutput:
        handle_de.close()
        handle_re.close()
        handle_sc.close()
        handle_dump.close()
        handle_single.close()
        handle_filt.close()


    # Write stats file
    handle_stat = open(outputDir + '/' + baseReadsFile + '.RSstat', 'w')
    handle_stat.write("## Hi-C processing\n")
    handle_stat.write("Valid_interaction_pairs\t" + str(valid_counter) + "\n")
    handle_stat.write("Valid_interaction_pairs_FF\t" + str(valid_counter_FF) + "\n")
    handle_stat.write("Valid_interaction_pairs_RR\t" + str(valid_counter_RR) + "\n")
    handle_stat.write("Valid_interaction_pairs_RF\t" + str(valid_counter_RF) + "\n")
    handle_stat.write("Valid_interaction_pairs_FR\t" + str(valid_counter_FR) + "\n")
    handle_stat.write("Dangling_end_pairs\t" + str(de_counter) + "\n")
    handle_stat.write("Religation_pairs\t" + str(re_counter) + "\n")
    handle_stat.write("Self_Cycle_pairs\t" + str(sc_counter) + "\n")
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

    handle_stat.close()

    if samOut:
        samfile.close()
