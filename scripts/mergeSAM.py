#!/usr/bin/python

## HiC-Pro
## Copyright (c) 2015-2016 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the
## terms of the BSD-3 licence.
## See the LICENCE file for details


"""
Script to pair 2 SAM/BAM files into one PE BAM
"""

# On 03/05/16 Ferhat made changes starting from
# ~/bin/HiC-Pro_2.7.2b/scripts/mergeSAM.py to make singletons possible to be
# reported

import getopt
import sys
import re
import pysam
from itertools import izip


def usage():
    """Usage function"""
    print "Usage : python mergeSAM.py"
    print "-f/--forward <forward read mapped file>"
    print "-r/--reverse <reverse read mapped file>"
    print "[-o/--output] <Output file. Default is stdin>"
    print "[-s/--single] <report singleton>"
    print "[-m/--multi] <report multiple hits>"
    print "[-q/--qual] <minimum reads mapping quality>"
    print "[-t/--stat] <generate a stat file>"
    print "[-v/--verbose] <Verbose>"
    print "[-h/--help] <Help>"
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "f:r:o:q:smtvh",
            ["forward=",
             "reverse=",
             "output=", "qual=",
             "single", "multi", "stat", "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts


def is_unique_bowtie2(read):
    ret = False
    if not read.is_unmapped and read.has_tag('AS'):
        if read.has_tag('XS'):
            primary = read.get_tag('AS')
            secondary = read.get_tag('XS')
            if (primary > secondary):
                ret = True
        else:
            ret = True

    return ret


# Remove everything after "/" in read's name
def get_read_name(read):
    name = read.qname
    return name.split("/", 1)[0]


def sam_flag(read1, read2, hr1, hr2):

    f1 = read1.flag
    f2 = read2.flag

    if not r1.is_unmapped:
        r1_chrom = hr1.getrname(r1.tid)
    else:
        r1_chrom = "*"
    if not r2.is_unmapped:
        r2_chrom = hr2.getrname(r2.tid)
    else:
        r2_chrom = "*"

    # Relevant bitwise flags (flag in an 11-bit binary number)
    # 1 The read is one of a pair
    # 2 The alignment is one end of a proper paired-end alignment
    # 4 The read has no reported alignments
    # 8 The read is one of a pair and has no reported alignments
    # 16 The alignment is to the reverse reference strand
    # 32 The other mate in the paired-end alignment is aligned
    #    to the reverse reference strand
    # 64 The read is the first (#1) mate in a pair
    # 128 The read is the second (#2) mate in a pair

    # The reads were mapped as single-end data, so should expect flags of 0
    # (map to the '+' strand) or 16 (map to the '-' strand) Output example: a
    # paired-end read that aligns to the reverse strand and is the first mate
    # in the pair will have flag 83 (= 64 + 16 + 2 + 1)

    if f1 & 0x4:
        f1 = f1 | 0x8

    if f2 & 0x4:
        f2 = f2 | 0x8

    if (not (f1 & 0x4) and not (f2 & 0x4)):
        # The flag should now indicate this is paired-end data
        f1 = f1 | 0x1
        f1 = f1 | 0x2
        f2 = f2 | 0x1
        f2 = f2 | 0x2

    # Indicate if the pair is on the reverse strand
    if f1 & 0x10:
        f2 = f2 | 0x20

    if f2 & 0x10:
        f1 = f1 | 0x20

    # Is this first or the second pair?
    f1 = f1 | 0x40
    f2 = f2 | 0x80

    # Insert the modified bitwise flags into the reads
    read1.flag = f1
    read2.flag = f2

    # Determine the RNEXT and PNEXT values (i.e. the positional values of a
    # read's pair) RNEXT
    if r1_chrom == r2_chrom:
        read1.rnext = r1.tid
        read2.rnext = r1.tid
    else:
        read1.rnext = r2.tid
        read2.rnext = r1.tid

    # PNEXT
    read1.pnext = read2.pos
    read2.pnext = read1.pos

    return(read1, read2)


if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    inputFile = None
    outputFile = None
    mapq = None
    report_single = False
    report_multi = False
    verbose = False
    stat = False
    output = "-"

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--forward"):
            R1file = arg
        elif opt in ("-r", "--reverse"):
            R2file = arg
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-q", "--qual"):
            mapq = arg
        elif opt in ("-s", "--single"):
            report_single = True
        elif opt in ("-m", "--multi"):
            report_multi = True
        elif opt in ("-t", "--stat"):
            stat = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print "## mergeBAM.py"
        print "## forward=", R1file
        print "## reverse=", R2file
        print "## output=", output
        print "## min mapq=", mapq
        print "## report_single=", report_single
        print "## report_multi=", report_multi
        print "## verbose=", verbose

    # Initialize variables
    tot_pairs_counter = 0
    multi_pairs_counter = 0
    uniq_pairs_counter = 0
    unmapped_pairs_counter = 0
    lowq_pairs_counter = 0
    multi_singles_counter = 0
    uniq_singles_counter = 0
    lowq_singles_counter = 0

    # local_counter = 0
    paired_reads_counter = 0
    singleton_counter = 0
    reads_counter = 0
    r1 = None
    r2 = None

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    if verbose:
        print "## Merging forward and reverse tags ..."

    with pysam.Samfile(R1file, "rb") as hr1, pysam.Samfile(R2file, "rb") as hr2:
        if output == "-":
            outfile = pysam.AlignmentFile(output, "w", template=hr1)
        else:
            outfile = pysam.AlignmentFile(output, "wb", template=hr1)
        for r1, r2 in izip(hr1.fetch(until_eof=True),
                           hr2.fetch(until_eof=True)):
            reads_counter += 1

            if (reads_counter % 1000000 == 0 and verbose):
                print "##", reads_counter

            if get_read_name(r1) == get_read_name(r2):
                # both unmapped
                if r1.is_unmapped and r2.is_unmapped:
                    unmapped_pairs_counter += 1
                    continue

                # both mapped
                elif not r1.is_unmapped and not r2.is_unmapped:
                    # quality
                    if (mapq is not None and
                        (r1.mapping_quality < int(mapq) or
                         r2.mapping_quality < int(mapq))):
                        lowq_pairs_counter += 1
                        continue

                    # Unique mapping
                    if is_unique_bowtie2(r1) and is_unique_bowtie2(r2):
                        uniq_pairs_counter += 1
                    else:
                        multi_pairs_counter += 1
                        if not report_multi:
                            continue
                # one end mapped, other is not
                else:
                    singleton_counter += 1
                    if not report_single:
                        continue
                    # first end is mapped, second is not
                    if not r1.is_unmapped:
                        # quality
                        if (mapq is not None) and \
                           (r1.mapping_quality < int(mapq)):
                            lowq_singles_counter += 1
                            continue

                        # Unique mapping
                        if is_unique_bowtie2(r1):
                            uniq_singles_counter += 1
                        else:
                            multi_singles_counter += 1
                            if not report_multi:
                                continue

                    else:  # second end is mapped, first is not
                        # quality
                        if (mapq is not None) and \
                           (r2.mapping_quality < int(mapq)):
                            lowq_singles_counter += 1
                            continue

                        # Unique mapping
                        if is_unique_bowtie2(r2):
                            uniq_singles_counter += 1
                        else:
                            multi_singles_counter += 1
                            if not report_multi:
                                continue

                tot_pairs_counter += 1
                (r1, r2) = sam_flag(r1, r2, hr1, hr2)

                # Write output
                outfile.write(r1)
                outfile.write(r2)

            else:
                print "Forward and reverse reads not paired."
                print "Check that BAM files are sorted."
                sys.exit(1)

    if stat:
        if output == '-':
            statfile = "pairing.stat"
        else:
            statfile = re.sub('.bam', '.pairstat', output)
        handle_stat = open(statfile, 'w')

        handle_stat.write(
            "Total_pairs_processed\t" +
            str(reads_counter) + "\t" +
            str(round(float(reads_counter)/reads_counter*100, .3)) +
            "\n")
        handle_stat.write(
            "Unmapped_pairs\t" +
            str(unmapped_pairs_counter) + "\t" +
            str(round(float(unmapped_pairs_counter)/reads_counter*100, 3)) +
            "\n")
        handle_stat.write(
            "Pairs_with_Singleton\t" +
            str(singleton_counter) + "\t" +
            str(round(float(singleton_counter)/reads_counter*100, 3)) + "\n")
        handle_stat.write(
            "Low_qual_pairs\t" + str(lowq_pairs_counter) + "\t" +
            str(round(float(lowq_pairs_counter)/reads_counter*100, 3)) + "\n")
        handle_stat.write(
            "Unique_paired_alignments\t" + str(uniq_pairs_counter) + "\t" +
            str(round(float(uniq_pairs_counter)/reads_counter*100, 3)) + "\n")
        handle_stat.write(
            "Multiple_pairs_alignments\t" + str(multi_pairs_counter) + "\t" +
            str(round(float(multi_pairs_counter)/reads_counter*100, 3)) + "\n")
        handle_stat.write(
            "Reported_pairs\t" + str(tot_pairs_counter) + "\t" +
            str(round(float(tot_pairs_counter)/reads_counter*100, 3)) + "\n")
        handle_stat.write(
            "Low_qual_singles\t" +
            str(lowq_singles_counter) + "\t" +
            str(round(float(lowq_singles_counter)/reads_counter*100, 3)) +
            "\n")
        handle_stat.write(
            "Unique_singles_alignments\t" +
            str(uniq_singles_counter) + "\t" +
            str(round(float(uniq_singles_counter)/reads_counter*100, 3)) +
            "\n")
        handle_stat.write(
            "Multiple_singles_alignments\t" +
            str(multi_singles_counter) + "\t" +
            str(round(float(multi_singles_counter)/reads_counter*100, 3)) +
            "\n")

        handle_stat.close()

    hr1.close()
    hr2.close()
    outfile.close()

