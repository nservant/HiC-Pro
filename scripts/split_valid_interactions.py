#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

"""
Script to split valid interactions into G1/G2 interaction files and calculate statistics
"""

import getopt
import sys
import os
import re
import pysam

def usage():
    """Usage function"""
    print("Usage : python split_valid_interactions.py")
    print("-i/--input <valid interaction file>")
    print("[-s/--stats] <stats file>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return


def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:s:vh",
            ["input=", "stats=",
             "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts


if __name__ == "__main__":
    ## Read command line arguments
    opts = get_args()
    inputfile = None
    statsFile = None
    verbose = False

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-s", "--stats"):
            statsFile = arg
        elif opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    ## Verbose mode
    if verbose:
        print("## split_valid_interactions.py")
        print("## input=", inputfile)
        print("## statsFile=", statsFile)
        print("## verbose=", verbose)

    ## AS counter
    vp_counter = 0
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

    G1cis_s = 0
    G1cis_l = 0
    G1trans = 0
    G2cis_s = 0
    G2cis_l = 0
    G2trans = 0

    ## Init output
    handle_g1 = open(inputfile.replace(".allValidPairs", "_G1.allValidPairs"), 'w')
    handle_g2 = open(inputfile.replace(".allValidPairs", "_G2.allValidPairs"), 'w')

    if verbose:
        print("## Splitting valid pairs interactions ...")
  
    with open(inputfile) as hr:
        for line in hr:
            isG1 = False
            isG2 = False

            vp_counter += 1
            h = line.rstrip().split("\t")
            haplotype = h[len(h)-1].split("-") ## always last column

            r1as = int(haplotype[0])
            r2as = int(haplotype[1])
            chr1 = h[1]
            chr2 = h[4]

            ## counter
            if r1as == 1 and r2as == 1:
                isG1 = True
                G1G1_ascounter += 1
                handle_g1.write(line)
            elif r1as == 1 and r2as == 0:
                isG1 = True
                G1U_ascounter += 1
                handle_g1.write(line)
            elif r1as == 0 and r2as == 1:
                isG1=True
                UG1_ascounter += 1
                handle_g1.write(line)
            
            elif r1as == 2 and r2as == 2:
                isG2 = True
                G2G2_ascounter += 1
                handle_g2.write(line)
            elif r1as == 2 and r2as == 0:
                isG2 = True
                G2U_ascounter += 1
                handle_g2.write(line)
            elif r1as == 0 and r2as == 2:
                isG2 = True
                UG2_ascounter += 1
                handle_g2.write(line)
            
            elif r1as == 1 and r2as == 2:
                G1G2_ascounter += 1
            elif r1as == 2 and r2as == 1:
                G2G1_ascounter += 1
            
            elif r1as == 3 or r2as == 3:
                CF_ascounter += 1
            else:
                UU_ascounter += 1

            ##Sats on distance
            if isG1 == True:
                if chr1 == chr2:
                    d = abs(int(h[5]) - int(h[2]))
                    if d <= 20000:
                        G1cis_s += 1
                    else:
                        G1cis_l += 1
                else:
                    G1trans += 1
            elif isG2 == True:
                if chr1 == chr2:
                    d = abs(int(h[5]) - int(h[2]))
                    if d <= 20000:
                        G2cis_s += 1
                    else:
                        G2cis_l += 1
                else:
                    G2trans += 1
                
            if (vp_counter % 100000 == 0 and verbose):
                print("##", vp_counter)

    if statsFile is not None:
        handle_stat = open(statsFile, 'w')            
        handle_stat.write("## HiC-Pro\n")
        handle_stat.write("## Allele specific information\n")
        handle_stat.write("Valid_pairs\t" + str(vp_counter) + "\n")
        handle_stat.write("Valid_pairs_from_ref_genome_(1-1)\t" + str(G1G1_ascounter) + "\n")
        handle_stat.write("Valid_pairs_from_ref_genome_with_one_unassigned_mate_(0-1/1-0)\t" + str(UG1_ascounter+G1U_ascounter) + "\n")
        handle_stat.write("Valid_pairs_from_alt_genome_(2-2)\t" + str(G2G2_ascounter) + "\n")
        handle_stat.write("Valid_pairs_from_alt_genome_with_one_unassigned_mate_(0-2/2-0)\t" + str(UG2_ascounter+G2U_ascounter) + "\n")
        handle_stat.write("Valid_pairs_from_alt_and_ref_genome_(1-2/2-1)\t" + str(G1G2_ascounter+G2G1_ascounter) + "\n")
        handle_stat.write("Valid_pairs_with_both_unassigned_mated_(0-0)\t" + str(UU_ascounter) + "\n")
        handle_stat.write("Valid_pairs_with_at_least_one_conflicting_mate_(3-)\t" + str(CF_ascounter) + "\n")
        handle_stat.write("cis_short_G1\t" + str(G1cis_s) + "\n")
        handle_stat.write("cis_long_G1\t" + str(G1cis_l) + "\n")
        handle_stat.write("trans_G1\t" + str(G1trans) + "\n")
        handle_stat.write("cis_short_G2\t" + str(G2cis_s) + "\n")
        handle_stat.write("cis_long_G2\t" + str(G2cis_l) + "\n")
        handle_stat.write("trans_G2\t" + str(G2trans) + "\n")
        handle_stat.close()



                   


     
