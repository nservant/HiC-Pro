#!/usr/bin/python

import getopt, sys, os
import re
from Bio import SeqIO

## Usage function
def usage():
    print "Usage : python readsTrim.py"
    print "-f/--fastq <Fastq file>"
    print "-c/--cutSite <Restriction enzyme cutting site>"
    print "[-v/--verbose <verbose mode>]"
    return


## getArgs function
def getArgs():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:c:vh", ["fastqFile=", "cutSite=", "verbose", "help"])
    except getopt.GetoptError as err:
	    usage()
	    sys.exit(-1)
    return opts


def trim_cutsite(fastq, cutsite):

    len_cs = len(cutsite) #cache this for later
    for record in fastq:
        index = record.seq.find(cutsite)
        if index == -1:
            #cutSite not found
            yield record
        else:
            #trim off the adaptor
            ##yield record[:index+len_cs]
            yield record[:index]



################################################################################
##
##  __MAIN__
##
################################################################################

## Read command line arguments
opts=getArgs()
fastqFile = None
cutSite = None
verbose = False

if len(opts)==0:
	usage()
	sys.exit()

for opt, arg in opts:
	if opt in ("-h", "--help"):
		usage()
		sys.exit()
	elif opt in ("-f", "--fastqFile"):
		fastqFile = arg
	elif opt in ("-c", "--cutSite"):
		cutSite = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
	else:
		assert False, "unhandled option"

## Verbose mode
if verbose:
    print "## readsTrim.py"
    print "## fastqFile=", fastqFile
    print "## cutSite=", cutSite
    print "## verbose=", verbose, "\n"
 
## Initialize
trim_counter = 0
total_counter = 0
trimmedFastq = re.sub(r'.fastq', '_trimmed.fastq', fastqFile)

## Read the fastq file
reads = SeqIO.parse(fastqFile, "fastq")
trimmed_reads = trim_cutsite(reads, cutSite)
count = SeqIO.write(trimmed_reads, trimmedFastq, "fastq") 
#print("Saved %i reads" % count)



