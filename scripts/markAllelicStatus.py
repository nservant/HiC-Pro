#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

"""
Script to assign an allelic status to a 'N' aligned BAM file
Note that the VCF file is loaded in RAM.
"""

import getopt
import sys
import os
import re
import pysam
import gzip

def usage():
    """Usage function"""
    print("Usage : python markAllelicStatus.py")
    print("-i/--ibam <BAM/SAM file of mapped reads>")
    print("-s/--snp <SNP file information - VCF format>")
    print("[-r/--rstat] <Generate a report with descriptive statistics>")
    print("[-t/--tag] <tag used to report allelic expression. Default XA>")
    print("[-o/--out] <output BAM file. Default is stdin>")
    print("[-v/--verbose] <verbose>")
    print("[-h/--help] <Help>")
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:s:t:o:rvh",
            ["ibam=",
             "snp=",
             "tag=",
             "output=",
             "rstat", "verbose", "help"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(-1)
    return opts


def get_snp_gt(gt, ref, alt): #!gt = ['0/1']         
    gtsnp = []
    
    ## gtsnp.append(ref)
    snp_geno = re.split('/|\|', gt)
    ## '.' are not considered
    if len(snp_geno) != 2:
        return [None, None]
    if (snp_geno[0] == '.') or (snp_geno[1] == '.'):
        return [None, None]

    ## First Allele
    if int(snp_geno[0]) == 0:
        gtsnp.append(ref)
    elif int(snp_geno[0]) == 1:
        gtsnp.append(alt)
    else:
        gtsnp.append(None)
    
    ## Second Allele
    if int(snp_geno[1]) == 0:
        gtsnp.append(ref)
    elif int(snp_geno[1]) == 1:
        gtsnp.append(alt)
    else:
        gtsnp.append(None)

    return gtsnp
    


def load_vcf(in_file, filter_qual=False, verbose=False, debug=False):
    """
    Load a VCF file in a dict object
    
    in_file = path to VCF file [character]
    filter_qual = if True, only good quality SNV are selected [boolean]
    verbose = if True, verbose mode[boolean]
    debug = if True, debug mode [boolean]
    """
    if verbose:
        print("## Loading VCF file {} ...".format(in_file))

    if in_file.endswith(".gz"):
        vcf_handle = gzip.open(in_file, 'rt')
    else:
        vcf_handle = open(in_file)    
    header = []
    samples = []
    snps = {}
    var_counter = 0
    snp_counter = 0
    for line in vcf_handle:
        line = line.rstrip()
        ## for now we don't care about the header
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            header = line.split('\t')
            header[0] = header[0][1:]
            samples = [s.split('.')[0] for s in header[9:]]
            if len(samples) > 1:
                print("Warning : Multisamples VCF detected. Only the first genotype will be used !", file=sys.stderr)
            continue
        else:
            fields = line.split('\t', 9)
            var_counter += 1
            n = len(fields)
            chrom, start, ref, alt, qfilter = fields[0], int(fields[1]) - 1, fields[3], fields[4], fields[6]
            ## Check format for first variant
            if var_counter == 1:
                formate = fields[8] if n>8 else None
                if formate.split(':')[0] != "GT":
                    print("Error : Invalid format - GT not detected at first position in ", file=sys.stderr)        
                    sys.exit(-1)

            genotypes  = fields[9].split('\t') if fields[9] else []
            geno = get_snp_gt(genotypes[0].split(':')[0], ref, alt)      
            if (filter_qual == False) or (filter_qual == True and qfilter == "PASS"):
                if debug:
                    print("{} - {} - {} -REF= {} -ALT= {} -G1= {} -G2= {}".format(str(chrom), str(start),
                          str(qfilter), str(ref), str(alt), str(geno[0]), str(geno[1])), file=sys.stderr)
                ## store only discriminant SNP
                if geno[0] != geno[1]:
                    snp_counter += 1
                    chrn = re.sub("^[Cc]hr","", chrom)
                    snps[(str(chrn), int(start), '1')] = geno[0]
                    snps[(str(chrn), int(start), '2')] = geno[1]

        if (var_counter % 100000 == 0 and verbose):
                print("##", var_counter)
 
    vcf_handle.close()
    if verbose:
        print("## Number of loaded SNPs = {} over {}".format(len(snps)/2, var_counter))
    return snps


def get_read_tag(read, tag):
    """
    Return a read's tag

    read  = read object [class pysam.AlignedSegment]      
    tag = name of the tag to return [character]
    """
    if t[0] == tag:
            return tag[1]
    return None


def get_mismatches_positions(read, base=None):
    """
    Extract the mismatch positions within the read for all mismatches or for a specified base
    Insertion are take into account and extracted from the CIGAR string
    Mismatch positions are extracted from the MD tag
    
    read = read object [class pysam.AlignedSegment]
    base = optional - base to look at [character]

    """
    md = read.get_tag('MD')
    x = -1 ## 0-based
    npos = [] ## 0-based
    digits = []    
    
    ## Get N pos in the read according to MD tag
    # for y in range(len(md)):
    y = 0
    while y < len(md):
        if md[y].isdigit():
            digits.append(md[y])
            #print "digits="+str(md[y])
        elif not md[y].isalnum(): ## simply ignore deletion
            y += 1
            while md[y].isalpha():
                #print "isAlpha="+md[y]
                y += 1
            #x-=1
            x += int(''.join(digits))
            digits = []
            continue
        elif md[y].isalpha():
            #print "alpha="+md[y]
            if len(digits) > 0:
                offset = int(''.join(digits))
                if base is None or (base is not None and md[y] == base):
                    npos.append(x + offset + 1)               
                digits = []
                x += offset + 1 
                #print "x pos="+str(x)
        y += 1

    #print npos
    ## Update N position if an insertion is detected upstream the N position
    ## l is the read based position
    if read.cigarstring.find("I") != -1:
        cig = read.cigartuples
        l = -1
        for t in cig:
            if t[0] == 1:
                for n in range(len(npos)):
                    #print "N=" + str(npos[n]) + " l=" + str(l) + "t=" + str(t)
                    if npos[n] > l:
                        npos[n] = npos[n] + t[1]
            if int(t[0]) != 3 and int(t[0]) != 2: ## skip splice junction
                l += t[1]
    return npos


def getGenomePos(read, pos):
    """
    Go from read position to genomic position

    read = read object [class pysam.AlignedSegment]
    pos = positions to convert [list]

    """

    ## Get genomic position, including Insertion/Deletion
    ngenomepos = []
    if len(pos) > 0:
        genomePos = read.get_reference_positions(full_length=True)
        for y in pos:
            if genomePos[y] == None:
                print("Warning : no genomic position found for {} at position {}".format(read.qname, y), file=sys.stderr)
                ngenomepos.append(None)
            else:
                ngenomepos.append(genomePos[y])
    return ngenomepos
    

def getBaseAt(read, pos):
    """
    Extract nucleotide within a read at a given set of positions

    read = read object [class pysam.AlignedSegment]
    pos = positions to convert [list]

    """
    nuc = [read.seq[p] for p in pos]
    return nuc


def getAllelicStatus(chrom, gpos, genotype, snps, debug=False):
    """
    For a given set of genomic position and associated genotype, compare to a snp file and return a code status
    0 : unassigned - no snp information extracted from the read
    1 : genotype from REF genome is found
    2 : genotype from ALT genome is found
    3 : conflicting information extracted from the read

    gpos = genomic position to look at [list]
    genotype = read genotype [list]
    snp = snp information extracted from VCF file
    """

    code = None
    g1_count, g2_count = 0, 0
    l = len(genotype)
    chrn = re.sub("^[Cc]hr","",chrom)

    for i in range(len(genotype)):
        #print >> sys.stderr, chrn, gpos[i], genotype[i]
        if gpos[i] != None:
            if ((str(chrn), int(gpos[i]), '1') in snps) and ((str(chrn), int(gpos[i]), '2') in snps):
                if snps[(str(chrn), int(gpos[i]), '1')] == genotype[i]:
                    g1_count += 1
                elif snps[(str(chrn), int(gpos[i]), '2')] == genotype[i]:
                    g2_count += 1
                else:
                    print("Warning : no SNPs found at position {} : {}. N ignored".format(chrom, str(gpos[i]+1)), file=sys.stderr)

    if g1_count > 0 and g2_count > 0:
        code = 3
    elif g1_count > 0 and g2_count == 0:
        code = 1
    elif g2_count > 0 and g1_count == 0:
        code = 2
    elif g1_count == 0 and g2_count == 0:
        code = 0

    return code

if __name__ == "__main__":

    # Read command line arguments
    opts = get_args()
    inputFile = None
    verbose = False
    report = False
    debug = False
    output = "-"
    tag = "XA"
    snps = {}
    
    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-s", "--snp"):
            snpFile = arg
        elif opt in ("-i", "--ibam"):
            mappedReadsFile = arg
        elif opt in ("-t", "--tag"):
            tag = arg
        elif opt in ("-o", "--out"):
            output = arg
        elif opt in ("-r", "--rstat"):
            report = True
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Initialize variables
    reads_counter = 0
    ua_counter = 0
    g1_counter = 0
    g2_counter = 0
    cf_counter = 0
    N_counter = 0

    # Read the SNP file
    snps = load_vcf(snpFile, filter_qual=False, verbose=verbose, debug=False)
    
    # Read the SAM/BAM file
    if verbose:
        print("## Opening SAM/BAM file {} ...".format(mappedReadsFile))
    infile = pysam.Samfile(mappedReadsFile, "rb")

    #samOut:
    if output == "-":
        outfile = pysam.AlignmentFile(output, "w", template=infile)
    else:
        outfile = pysam.AlignmentFile(output, "wb", template=infile)
     
   # Verbose mode                                                                                                                                                        
    if verbose:
        print("## " + __file__)
        print("## ibam=", mappedReadsFile)
        print("## snpFile=", snpFile)
        print("## tag=", tag)
        print("## output=" + output)
        print("## verbose={}\n".format(verbose))

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    if verbose:
        print("## Assigning allele specific information ...")

  
    for read in infile.fetch(until_eof=True):
        reads_counter += 1
        if not read.is_unmapped:## and read.cigarstring.find("D") != -1:
            read_chrom = infile.get_reference_name(read.tid)
            Nreadpos = get_mismatches_positions(read, base="N")
            if (len(Nreadpos) > 0):
                N_counter += len(Nreadpos)
                Ngenomepos = getGenomePos(read, Nreadpos)
                Nbase = getBaseAt(read, Nreadpos)
                tagval = getAllelicStatus(read_chrom, Ngenomepos, Nbase, snps, debug=debug)
                read.set_tag(tag, tagval)
            
                if debug:
                    for i in range(len(Nreadpos)):
                        if Ngenomepos[i] != None:
                            print("{}\t{}\t{}\t{}/N/{}\t{}".format(str(read_chrom), str(Ngenomepos[i]),
                                  str(Ngenomepos[i]+1), str(read.query_name), str(Nbase[i]), str(tagval)), file=sys.stderr)
                if tagval == 0:
                    ua_counter += 1
                elif tagval == 1:
                    g1_counter += 1
                elif tagval == 2:
                    g2_counter += 1
                elif tagval == 3:
                    cf_counter += 1                                        
            else:
                read.set_tag(tag, 0)
                ua_counter += 1
        else:
            read.set_tag(tag, 0)
            ua_counter += 1

        outfile.write(read)

        if (reads_counter % 100000 == 0 and verbose):
            print("##", reads_counter)

    # Close handler
 
    # Write stats file
    if report:
        handle_stat = open(re.sub(r'\.bam$|\.sam$', '.allelstat', output), 'w')
        handle_stat.write("## " + __file__ + "\n")
        handle_stat.write("## ibam=" + mappedReadsFile + "\n")
        handle_stat.write("## snpFile=" + snpFile + "\n")
        handle_stat.write("## tag=" + tag + "\n")
        handle_stat.write("## output=" + output + "\n")
        handle_stat.write("## verbose=" + str(verbose) + "\n")
        handle_stat.write("## =========================\n")

        handle_stat.write("Total number of snps loaded\t" + str(len(snps)/2) + "\n")
               
        handle_stat.write("## =========================\n")

        handle_stat.write("Total number of reads\t" + str(reads_counter) + "\t100" + "\n")
        handle_stat.write("Number of reads with at least one 'N'\t" + str(N_counter) + "\t" + str(round(float(N_counter)/int(reads_counter)*100,3)) + "\n")
        handle_stat.write("Number of reads assigned to ref genome\t" + str(g1_counter) + "\t" + str(round(float(g1_counter)/int(reads_counter)*100,3)) + "\n")
        handle_stat.write("Number of reads assigned to alt genome\t" + str(g2_counter) + "\t" + str(round(float(g2_counter)/int(reads_counter)*100,3)) + "\n")
        handle_stat.write("Number of conflicting reads\t" + str(cf_counter) + "\t" + str(round(float(cf_counter)/int(reads_counter)*100,3)) + "\n")
        handle_stat.write("Number of unassigned reads\t" + str(ua_counter) + "\t" + str(round(float(ua_counter)/int(reads_counter)*100,3)) + "\n")
        handle_stat.close()
                
    infile.close()
    outfile.close()
