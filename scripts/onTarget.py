#!/usr/bin/env python

# HiC-Pro
# Copyleft 2016 Institut Curie
# Author(s): Nicolas Servant, Eric Viara
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the BSD License

"""
Script to keep valid 3C products overlapping with Target
The list of valid pairs on targets are printed on stdout
Errors and verbose in stderr
"""
import getopt
import sys
from bx.intervals.intersection import Intersecter, Interval

def usage():
    """Usage function"""
    print("Usage : python onTarget.py")
    print("-i/--inFile <Valid Pairs file>")
    print("-t/--target <BED file of targets>")
    print("[-s/--stats] <Stats file>")
    print("[-c/--cis] <Report only capture-capture interactions. Otherwise both capture-capture and capture-reporter interactions are returned>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:t:s:cvh",
            ["inFile=",
             "target=",
             "stats=",
             "cis", "verbose"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

def load_bed(in_file, verbose=False):
    """
    Read a BED file and store the intervals in a tree

    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    in_file = input file [character]
    verbose = verbose mode [logical]

    """
    intervals = {}
    if verbose:
        print("## Loading BED file {} ...".format(in_file), file=sys.stderr)
    with open(in_file, 'r') as bed_handle:
        nline = 0
        for line in bed_handle:
            nline +=1
            bedtab = line.strip().split("\t")
            try:
                chromosome, start, end = bedtab[:3]
            except ValueError:
                print("Warning : wrong input format in line {}. Not a BED file !?".format(nline),
                      file=sys.stderr)
                continue

            # BED files are zero-based as Intervals objects
            start = int(start)  # + 1
            end = int(end)
            fragl = abs(end - start)
        
            if chromosome in intervals:
                tree = intervals[chromosome]
                tree.add_interval(Interval(start, end))
            else:
                tree = Intersecter()
                tree.add_interval(Interval(start, end))
                intervals[chromosome] = tree
    
    return intervals

if __name__ == "__main__":
    # Read command line arguments
    opts = get_args()
    verbose = False
    cis = False
    statsFile = None

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--inFile"):
            inFile = arg
        elif opt in ("-t", "--target"):
            target = arg
        elif opt in ("-c", "--cis"):
            cis = True
        elif opt in ("-s", "--stats"):
            statsFile = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # Verbose mode
    if verbose:
        print("## onTarget.py", file=sys.stderr)
        print("## inFile ={}".format(inFile), file=sys.stderr)
        print("## target ={}".format(target), file=sys.stderr)
        print("## cis ={}".format(cis), file=sys.stderr)
        print("## verbose ={}\n".format(verbose), file=sys.stderr)

    # Initialize variables
    vp_counter = 0
    ontarget_counter = 0
    ontarget_cap_cap_counter = 0
    ontarget_cap_rep_counter = 0

    # Read the BED file
    targetInter = load_bed(target, verbose)
        
    # Read the SAM/BAM file
    if verbose:
        print("## Opening valid pairs file {}...".format(inFile), file=sys.stderr)
    
    vp_handle = open(inFile)

    with open(inFile) as vp_handle:
        for line in vp_handle:
            vp_counter += 1
            sline = line.strip().split("\t")
            try:
                chr1 = sline[1]
                pos1 = int(sline[2])
                chr2 = sline[4]
                pos2 = int(sline[5])
            except ValueError:
                print("Warning : wrong input format in line {}. Not a BED file !?".format(nline), file=sys.stderr)
                continue
            
            res1 = []
            res2 = []
            if chr1 in targetInter:
                res1 = targetInter[chr1].find(pos1, pos1+1)
            if chr2 in targetInter:
                res2 = targetInter[chr2].find(pos2, pos2+1)
                
            if len(res1) > 0 and len(res2) > 0:
                ontarget_counter += 1
                ontarget_cap_cap_counter += 1
                print(line.strip())
            elif len(res1) > 0 or len(res2) > 0:
                ontarget_counter += 1
                ontarget_cap_rep_counter += 1
                if not cis:
                    print(line.strip())

        if statsFile is not None:
            if verbose:
                print("## Writing stats in {}...".format(statsFile), file=sys.stderr)
            f = open(statsFile, 'a')
            f.write("valid_pairs_on_target\t" + str(ontarget_counter) + "\n")
            f.write("valid_pairs_on_target_cap_cap\t" + str(ontarget_cap_cap_counter) + "\n")
            f.write("valid_pairs_on_target_cap_rep\t" + str(ontarget_cap_rep_counter) + "\n")
            f.closed

