#!/usr/bin/python

import getopt, sys, os
import re
import pysam
from bx.intervals.intersection import Intersecter, Interval

## Usage function
def usage():
    print "Usage : python overlapMapped2HiCFragments.py"
    print "-f/--fragmentFile <Restriction fragment file GFF3>"
    print "-r/--mappedReadsFile <BAM/SAM file of mapped reads>"
    print "[-o/--outputDir] <Output directory. Default is current directory>"
    print "[-s/--shortestInsertSize] <Shortest insert size of mapped reads to consider>"
    print "[-l/--longestInsertSize] <Longest insert size of mapped reads to consider>"
    print "[-a/--all] <Write all additional output files, with information about the discarded reads>"
    print "[-v/--verbose] <Verbose>"
    print "[-h/--help] <Help>"
    return


## getArgs function
def getArgs():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:r:o:s:l:vah", ["fragmentFile=", "mappedReadsFile=", "outputDir=", "minInsertSize=", "maxInsertSize", "verbose", "all", "help"])
    except getopt.GetoptError as err:
	    usage()
	    sys.exit(-1)
    return opts


## getReadStrand
## Conversion of read position to naive strand representation
##
## read = [AlignedRead]
def getReadStrand ( read ):
    strand = "+"
    if read.is_reverse:
        strand="-"
    return strand


## getReadStart
## Return read 5' end (zero-based) according to the mapping strand
##
## read = [AlignedRead]
def getRead5end ( read ):
    if (read.is_reverse):
        pos = read.pos + read.alen + 1 ##(50 - 5) + 1 ## zero-based transformation
    else:
        pos = read.pos 
    return pos

## getOrderedReads
## The sequecing is usually not oriented. Reorient the reads so that r1 is always before r2
##
## read1 = [AlignedRead]
## read2 = [AlignedRead]
def getOrderedReads ( read1, read2 ):
    if read1.tid == read2.tid:
        if getRead5end(read1) < getRead5end(read2):
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
  
## loadRestrictionFragment
## Read a BED file and store the intervals in a tree
## Intervals are zero-based objects. The output object is a hash table with one search tree per chromosome
##
## in_file = input file [character]
## verbose = verbose mode [logical]
def loadRestrictionFragment( in_file, verbose ):
    resFrag={}
    if verbose:
        print "## Loading Restriction File Intervals '", in_file,"'..."

    bed_handle = open(in_file)
    for line in bed_handle:
        bedtab = line.split("\t")
        number_of_columns = len(bedtab)
        try:
            chromosome, start, end, name = bedtab[:4]
        except ValueError:
            print error_message.format(count+1, number_of_columns, line.strip())
            continue

        ## BED files are zero-based as Intervals objects
        start = int(start) ## + 1
        end=int(end)
        name=name.strip()
        if chromosome in resFrag.keys():
            tree = resFrag[chromosome]
            tree.add_interval( Interval(start, end, value={'name':name}) )
        else:
            tree = Intersecter()
            tree.add_interval( Interval(start, end, value={'name':name}) )
            resFrag[chromosome]=tree
    bed_handle.close()
    return resFrag


## getOverlappingRestrictionFragment
## Intersect a given read with the set of restriction fragments
##
## resFrag = the restriction fragments [hash] 
## chrom = the chromosome to look at [character] 
## read = the read to intersect [AlignedRead]
def getOverlappingRestrictionFragment ( resFrag, chrom, read ):
    ## Get 5' end
    pos = getRead5end(read)
    ## Overlap with the 5' end of the read (zero-based)
    resfrag = resFrag[chrom].find(pos, pos+1)
    if len(resfrag)>1:
        print "Error : ",len(resfrag)," restriction fragments found for ",read.qname
        
    return resfrag[0]


## isSelfCircle
## Both reads are expected to be on the same restriction fragments
## Check the orientation of reads <-->
## read1 : [AlignedRead]
## read2 : [AlignedRead]
def isSelfCircle ( read1, read2 ):
    ret = False;
    ## Get oriented reads
    r1, r2 = getOrderedReads(read1, read2)
    ## 1<- ->2 or 2<- ->1  
    if getReadStrand(r1) == "-" and getReadStrand(r2) == "+":
        ret=True
    return ret


## isDanglingEnd
## Both reads are expected to be on the same restriction fragments
## Check the orientation of reads -><-
## read1 : [AlignedRead]
## read2 : [AlignedRead]
def isDanglingEnd ( read1, read2 ):
    ret = False;
    ## Get oriented reads
    r1, r2 = getOrderedReads(read1, read2)
    ## 1-> <-2 or 2-> <-1
    if getReadStrand(r1) == "+" and getReadStrand(r2) == "-":
        ret=True
    return ret
   

## getValidOrientation
## Both reads are expected to be on the different restriction fragments
## Check the orientation of reads ->-> / <-<- / -><- / <-->
## read1 : [AlignedRead]
## read2 : [AlignedRead]
def getValidOrientation ( read1, read2 ):    
    ## Get oriented reads
    r1, r2 = getOrderedReads(read1, read2)

    direction=None
    if getReadStrand(r1) == "+" and getReadStrand(r2) == "+":
        direction="FF"
    elif  getReadStrand(r1) == "-" and getReadStrand(r2) == "-":
        direction="RR"
    elif  getReadStrand(r1) == "+" and getReadStrand(r2) == "-":
        direction="FR"
    elif  getReadStrand(r1) == "-" and getReadStrand(r2) == "+":
        direction="RF"
    
    return direction


## getPEFragmentSize
## Calculte the size of the DNA fragment library
## read1 : [AlignedRead]
## read2 : [AlignedRead]    
## resfrag1 = restrictin fragment overlapping the R1 read [interval]
## resfrag1 = restrictin fragment overlapping the R1 read [interval]
## interactionType = Type of interaction from getInteractionType() [character]
def getPEFragmentSize ( read1, read2, resFrag1, resFrag2, interactionType ):
    fragmentsize=None
    ## Get oriented reads
    r1, r2 = getOrderedReads(read1, read2)
    if r1==read2:
        rfrag1=resFrag2
        rfrag2=resFrag1
    else:
        rfrag1=resFrag1
        rfrag2=resFrag2

    r1pos = getRead5end(r1)
    r2pos = getRead5end(r2)

    if interactionType == "DE":
        fragmentsize=r2pos - r1pos 
    elif interactionType == "SC":
        fragmentsize=(r1pos - rfrag1.start) + (rfrag2.end - r2pos)
    elif interactionType == "VALID":
        if getReadStrand(r1)=="+":
            dr1 = rfrag1.end - r1pos
        else:
            dr1 = r1pos - rfrag1.start
        if getReadStrand(r2)=="+":
            dr2 = rfrag2.end - r2pos
        else:
            dr2 = r2pos - rfrag2.start
        fragmentsize = dr2 + dr1

    return fragmentsize


## getInteractionType
## For a given reads pair and their related restriction fragment, classify the 3C products as :
## - Interaction
## - Self circle
## - Dangling end
## - Unknown
##
## read1 = the R1 read of the pair [AlignedRead]
## read1_chrom = the chromosome of R1 read [character]
## resfrag1 = restrictin fragment overlapping the R1 read [interval]
## read2 = the R2 read of the pair [AlignedRead]
## read2_chrom = the chromosome of R2 read [character]
## resfrag2 = restrictin fragment overlapping the R2 read [interval]
## verbose = verbose mode [logical]
def  getInteractionType( read1, read1_chrom, resfrag1, read2, read2_chrom, resfrag2, verbose):               
    ## If returned InteractionType=None -> Same restriction fragment and same strand = Error
    InteractionType = None
    strand1 = getReadStrand(read1)
    strand2 = getReadStrand(read2)
    chr1 = read1_chrom
    chr2 = read2_chrom
        
    ##print chr1,"\t",read1.pos+1,"\t",strand1,"\t",resfrag1.value['name'],"|mm9|",chr1,":",resfrag1.start+1,"-",resfrag1.end,"\t",chr2,"\t",read2.pos+1,"\t",strand2,"\t",resfrag2.value['name'],"|mm9|",chr2,":",resfrag2.start+1,"-",resfrag2.end,
        
    if not (r1.is_unmapped) and  not (r2.is_unmapped):
        ## same restriction fragment
        if resfrag1 == resfrag2:
            ## Self_circle <- ->
            if isSelfCircle(read1, read2):
                InteractionType="SC"
            ## Dangling_end -> <-
            elif isDanglingEnd(read1, read2):
                InteractionType="DE"
        else:
            InteractionType="VALID"
    elif r1.is_unmapped or r2.is_unmapped:
        InteractionType="SINGLE"
                    
    return InteractionType




################################################################################
##
##  __MAIN__
##
################################################################################

## Read command line arguments
opts=getArgs()
inputFile = None
outputFile = None
verbose = False
allOutput = False
minInsertSize = None
maxInsertSize = None
outputDir = "."

if len(opts)==0:
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
	elif opt in ("-v", "--verbose"):
		verbose = True	
	else:
		assert False, "unhandled option"

## Verbose mode
if verbose:
    print "## overlapMapped2HiCFragments.py"
    print "## mappedReadsFile=", mappedReadsFile
    print "## fragmentFile=", fragmentFile
    print "## minInsertSize=", minInsertSize
    print "## maxInsertSize=", maxInsertSize
    print "## allOuput=", allOutput
    print "## verbose=", verbose, "\n"
 
## Initialize variables
reads_counter=0
de_counter=0
sc_counter=0
valid_counter=0
valid_counter_FF=0
valid_counter_RR=0
valid_counter_FR=0
valid_counter_RF=0
single_counter=0
error_counter=0
##dist_hash={}

baseReadsFile = os.path.basename(mappedReadsFile)
baseReadsFile = re.sub(r'.bam|.sam', '', baseReadsFile)

## Open handlers for output files
handle_valid = open(outputDir + '/' + baseReadsFile + '.validPairs', 'w')

if allOutput:
    handle_de = open(outputDir + '/' + baseReadsFile + '.DEPairs', 'w')
    handle_sc = open(outputDir + '/' + baseReadsFile + '.SCPairs', 'w')
    handle_err = open(outputDir + '/' + baseReadsFile + '.ErrorPairs', 'w')
    handle_single = open(outputDir + '/' + baseReadsFile + '.SinglePairs', 'w')

## Read the BED file
resFrag=loadRestrictionFragment(fragmentFile, verbose)

## Read the SAM/BAM file
if verbose:
    print "## Opening SAM/BAM file '", mappedReadsFile,"'..."

samfile = pysam.Samfile( mappedReadsFile, "r" )

## Reads are 0-based too (for both SAM and BAM format)
## Loop on all reads
if verbose:
    print "## Classifying Interactions ..."

for read in samfile.fetch():
    reads_counter += 1
    cur_handler=None

    ## First mate
    if read.is_read1:
        r1 = read
        r1_chrom = samfile.getrname(r1.tid)
        r1_resfrag = getOverlappingRestrictionFragment(resFrag, r1_chrom, r1)
  
    ## Second mate
    elif read.is_read2:
        r2 = read
        r2_chrom = samfile.getrname(r2.tid)
        r2_resfrag = getOverlappingRestrictionFragment(resFrag, r2_chrom, r2)
        interactionType=getInteractionType(r1, r1_chrom, r1_resfrag, r2, r2_chrom, r2_resfrag, verbose)
        dist=getPEFragmentSize(r1, r2, r1_resfrag, r2_resfrag, interactionType)
        if (minInsertSize != None and dist != None and dist < int(minInsertSize)) or (maxInsertSize != None and dist != None and dist > int(maxInsertSize)):
            interactionType=None
        ##else:
        ##    dist_hash[str(dist)]=dist_hash[str(dist)]+1 if dist_hash.has_key(str(dist)) else 1

        if interactionType == "VALID":
            valid_counter+=1
            cur_handler = handle_valid
            validType=getValidOrientation(r1, r2)
            if validType == "RR":  valid_counter_RR+=1
            elif validType == "FF":  valid_counter_FF+=1
            elif validType == "FR":  valid_counter_FR+=1
            elif validType == "RF":  valid_counter_RF+=1

        elif interactionType == "DE":
            de_counter+=1
            cur_handler = handle_de if allOutput else None

        elif interactionType == "SC":
            sc_counter+=1
            cur_handler = handle_sc if allOutput else None
           
        elif interactionType == "SINGLE":
            single_counter+=1
            cur_handler = handle_single if allOutput else None

        else:
            error_counter+=1
            cur_handler = handle_err if allOutput else None

        if cur_handler != None:
            cur_handler.write(r1.qname + "\t" + r1_chrom + "\t" + str(r1.pos) + "\t" + str(getReadStrand(r1)) + "\t" + r2_chrom + "\t" + str(r2.pos) + "\t" + str(getReadStrand(r2)) + "\t" + str(dist) + "\n")

    if (reads_counter % 100000 == 0 and verbose):
        print "##", reads_counter
 
## Close handler
handle_valid.close()
if allOutput:
    handle_de.close()
    handle_sc.close()
    handle_err.close()
    handle_single.close()

if verbose:
    print "Valid interaction pairs\t", valid_counter
    print "Valid interaction pairs FF\t", valid_counter_FF
    print "Valid interaction pairs RR\t", valid_counter_RR
    print "Valid interaction pairs RF\t", valid_counter_RF
    print "Valid interaction pairs FR\t", valid_counter_FR
    print "Dangling end pairs\t", de_counter
    print "Self Cycle pairs\t", sc_counter
    print "Single-end pairs\t", single_counter
    print "Error pairs\t", error_counter, "\n"


samfile.close()
