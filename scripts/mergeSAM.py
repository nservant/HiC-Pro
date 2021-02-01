#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details


"""
Script to pair 2 SAM/BAM files into one PE BAM
- On 03/05/16 Ferhat made changes starting from ~/bin/HiC-Pro_2.7.2b/scripts/mergeSAM.py 
to make singletons possible to be reported
"""

import getopt
import sys
import os
import re
import pysam

def usage():
    """Usage function"""
    print("Usage : python mergeSAM.py")
    print("-f/--forward <forward read mapped file>")
    print("-r/--reverse <reverse read mapped file>")
    print("[-o/--output] <Output file. Default is stdin>")
    print("[-s/--single] <report singleton>")
    print("[-m/--multi] <report multiple hits>")
    print("[-q/--qual] <minimum reads mapping quality>")
    print("[-t/--stat] <generate a stat file>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
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
            primary =  read.get_tag('AS')
            secondary = read.get_tag('XS')
            if (primary > secondary):
                ret = True
        else:
            ret = True
    return ret

## Remove everything after "/" or " " in read's name
def get_read_name(read):
    name = read.query_name
    #return name.split("/",1)[0]
    return re.split('/| ', name)[0]

def sam_flag(read1, read2, hr1, hr2):
	
    f1 = read1.flag
    f2 = read2.flag

    if r1.is_unmapped == False:
        r1_chrom = hr1.get_reference_name(r1.reference_id)
    else:
        r1_chrom = "*"
    if r2.is_unmapped == False:
        r2_chrom = hr2.get_reference_name(r2.reference_id)
    else:
        r2_chrom="*"

    ##Relevant bitwise flags (flag in an 11-bit binary number)
    ##1 The read is one of a pair
    ##2 The alignment is one end of a proper paired-end alignment
    ##4 The read has no reported alignments
    ##8 The read is one of a pair and has no reported alignments
    ##16 The alignment is to the reverse reference strand
    ##32 The other mate in the paired-end alignment is aligned to the reverse reference strand
    ##64 The read is the first (#1) mate in a pair
    ##128 The read is the second (#2) mate in a pair
  
    ##The reads were mapped as single-end data, so should expect flags of 
    ##0 (map to the '+' strand) or 16 (map to the '-' strand)
    ##Output example: a paired-end read that aligns to the reverse strand 
    ##and is the first mate in the pair will have flag 83 (= 64 + 16 + 2 + 1)
  
    if f1 & 0x4:
        f1 = f1 | 0x8

    if f2 & 0x4:
        f2 = f2 | 0x8
    
    if (not (f1 & 0x4) and not (f2 & 0x4)):
        ##The flag should now indicate this is paired-end data
        f1 = f1 | 0x1
        f1 = f1 | 0x2
        f2 = f2 | 0x1
        f2 = f2 | 0x2  
    
    ##Indicate if the pair is on the reverse strand
    if f1 & 0x10:
        f2 = f2 | 0x20
  
    if f2 & 0x10:
        f1 = f1 | 0x20
  
    ##Is this first or the second pair?
    f1 = f1 | 0x40
    f2 = f2 | 0x80
  
    ##Insert the modified bitwise flags into the reads
    read1.flag = f1
    read2.flag = f2
	
    ##Determine the RNEXT and PNEXT values (i.e. the positional values of a read's pair)
    #RNEXT
    if r1_chrom == r2_chrom:
        read1.next_reference_id = r1.reference_id
        read2.next_reference_id = r1.reference_id
    else:
        read1.next_reference_id = r2.reference_id
        read2.next_reference_id = r1.reference_id
    #PNEXT
    read1.next_reference_start = read2.reference_start
    read2.next_reference_start = read1.reference_start

    return(read1, read2)



if __name__ == "__main__":
    ## Read command line arguments
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

    ## Verbose mode
    if verbose:
        print("## mergeBAM.py")
        print("## forward=", R1file)
        print("## reverse=", R2file)
        print("## output=", output)
        print("## min mapq=", mapq)
        print("## report_single=", report_single)
        print("## report_multi=", report_multi)
        print("## verbose=", verbose)

    ## Initialize variables
    tot_pairs_counter = 0
    multi_pairs_counter = 0
    uniq_pairs_counter = 0
    unmapped_pairs_counter = 0 
    lowq_pairs_counter = 0
    multi_singles_counter = 0
    uniq_singles_counter = 0
    lowq_singles_counter = 0

    #local_counter = 0
    paired_reads_counter = 0
    singleton_counter = 0
    reads_counter = 0
    r1 = None
    r2 = None

    ## Reads are 0-based too (for both SAM and BAM format)
    ## Loop on all reads
    if verbose:
        print("## Merging forward and reverse tags ...")
    
    with pysam.Samfile(R1file, "rb") as hr1, pysam.Samfile(R2file, "rb") as hr2: 
        if output == "-":
            outfile = pysam.AlignmentFile(output, "w", template=hr1)
        else:
            outfile = pysam.AlignmentFile(output, "wb", template=hr1)
	
        for r1, r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
            reads_counter +=1
            if (reads_counter % 1000000 == 0 and verbose):
                print("##", reads_counter)
                
            if get_read_name(r1) == get_read_name(r2):
                ## both unmapped
                if r1.is_unmapped == True and r2.is_unmapped == True:
                    unmapped_pairs_counter += 1
                    continue
                    
                ## both mapped
                elif r1.is_unmapped == False and r2.is_unmapped == False:
                    ## quality
                    if mapq != None and (r1.mapping_quality < int(mapq) or r2.mapping_quality < int(mapq)):
                        lowq_pairs_counter += 1
                        continue
                 
                    ## Unique mapping
                    if is_unique_bowtie2(r1) == True and is_unique_bowtie2(r2) == True:
                        uniq_pairs_counter += 1
                    else:
                        multi_pairs_counter += 1
                        if report_multi == False:
                            continue

                ## One mate maped
                else:
                    singleton_counter += 1
                    if report_single == False:
                        continue
                    if r1.is_unmapped == False:  ## first end is mapped, second is not
                        ## quality
                        if mapq != None and (r1.mapping_quality < int(mapq)): 
                            lowq_singles_counter += 1
                            continue
                        ## Unique mapping
                        if is_unique_bowtie2(r1) == True:
                            uniq_singles_counter += 1
                        else:
                            multi_singles_counter += 1
                            if report_multi == False:
                                continue
                    else:  ## second end is mapped, first is not
                        ## quality
                        if mapq != None and (r2.mapping_quality < int(mapq)): 
                            lowq_singles_counter += 1
                            continue
                        ## Unique mapping
                        if is_unique_bowtie2(r2) == True:
                            uniq_singles_counter += 1
                        else:
                            multi_singles_counter += 1
                            if report_multi == False:
                                continue

                tot_pairs_counter += 1          
                (r1, r2) = sam_flag(r1,r2, hr1, hr2)

                ## Write output
                outfile.write(r1)
                outfile.write(r2)
                
            else:
                print("Forward and reverse reads not paired. Check that BAM files have the same read names and are sorted.")
                sys.exit(1)

        if stat:
            if output == '-':
                statfile = "pairing.stat"
            else:
                statfile = re.sub('\.bam$', '.pairstat', output)
            with open(statfile, 'w') as handle_stat:
                handle_stat.write("Total_pairs_processed\t" + str(reads_counter) + "\t" + str(round(float(reads_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Unmapped_pairs\t" + str(unmapped_pairs_counter) + "\t" + str(round(float(unmapped_pairs_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Low_qual_pairs\t" + str(lowq_pairs_counter) + "\t" + str(round(float(lowq_pairs_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Unique_paired_alignments\t" + str(uniq_pairs_counter) + "\t" + str(round(float(uniq_pairs_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Multiple_pairs_alignments\t" + str(multi_pairs_counter) + "\t" + str(round(float(multi_pairs_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Pairs_with_singleton\t" + str(singleton_counter) + "\t" + str(round(float(singleton_counter)/float(reads_counter)*100,3)) + "\n")  
                handle_stat.write("Low_qual_singleton\t" + str(lowq_singles_counter) + "\t" + str(round(float(lowq_singles_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Unique_singleton_alignments\t" + str(uniq_singles_counter) + "\t" + str(round(float(uniq_singles_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Multiple_singleton_alignments\t" + str(multi_singles_counter) + "\t" + str(round(float(multi_singles_counter)/float(reads_counter)*100,3)) + "\n")
                handle_stat.write("Reported_pairs\t" + str(tot_pairs_counter) + "\t" + str(round(float(tot_pairs_counter)/float(reads_counter)*100,3)) + "\n")
    hr1.close()
    hr2.close()
    outfile.close()

