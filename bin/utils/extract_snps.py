#!/usr/bin/env python
# Author(s): Nicolas Servant
# Contact: nicolas.servant@curie.fr
# This software is distributed without any guarantee under the terms of the GNU General 
# Public License, either Version 2, June 1991 or Version 3, June 2007.

"""
Script to extract informative SNPs from a multisamples VCF
If no --ref is specified, use the REF allele as is. Otherwise, replace the REF allele by the specified one.
Phased data are not considered

Example :
./extract_snps.py -i mgp.v2.snps.annot.reformat.vcf -r CASTEiJ -a FVB_NJ -q -v > mgp.v2.snps.annot.reformat_CAST_FVB_PASS.vcf
"""

import getopt
import sys
import os
import pysam
import re
import gzip


def usage():
    print("This script was designed to extract informative SNPs information from two parental genotypes, and return the F1 genotype.")

    """Usage function"""
    print("Usage : python extract_snps.py")
    print("-i/--vcf <input VCF file information>")
    print("-a/--alt <sample name for alternative allele>")
    print("[-r/--ref] <sample name for reference allele. Default is the same as in the initial file>")
    print("[-f/--filt] <Filtering level. 0 - no filtering. 1 - based on FI information. 2 - based on FILTER information. Default is 2>")
    print("[-x/--exclude] <Exclude potential contaminants from the list of SNPs which are equal to ALT genotype, i.e REF is expected to be contaminated>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:r:a:f:x:vh",
            ["vcf=",
             "alt=",
             "ref=",
             "filt=", 
             "excl=",
             "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

## 0/0, 1/1, A, T, 1
## 1/1, 2/2, A, T, C, 1
def get_filter_snp_gt(gref, galt, ref, alt, conta):

    #print gref
    #print galt
    #print ref
    #print alt

    ref_geno = re.split('/|\|', gref)
    alt_geno = re.split('/|\|', galt)
        
    ## '.' are not considered                                                                                                                                    
    if len(alt_geno) != 2 or len(ref_geno) != 2:
        return -1
    elif ref_geno[0] == "." or ref_geno[1] == "." or alt_geno[0] == "." or alt_geno[1] == ".":
        return -1
    ## remove heterogygous
    elif ref_geno[0] != ref_geno[1] or alt_geno[0] != alt_geno[1]:
        return -2
    else:
        ## Keep homozyguous SNPs
        ref_snp = ref_geno[0]
        alt_snp = alt_geno[0]

        ## Non informative SNPs
        if ref_snp == alt_snp:
            return -3
        ## Remove SNPs that are different from the reference and therefore likely to be wrongly assigned to the alternative
        elif len(conta)>0:
            for i in range(len(conta)):
                conta_geno = re.split('/|\|', conta[i])
                ## if conta = alt remove snps
                if conta_geno[0] == alt_geno[0] or conta_geno[1] == alt_geno[1]:
                    return -4
        
        ## check alt and ref alleles
        alleles = []
        alleles.append(ref)
        for a in alt.split(','):
            alleles.append(a)
            
        return [alleles[int(ref_snp)], alleles[int(alt_snp)]]


class Check_input_parameters:
    def ref_integrity(self, refSample, refidx):
        if refSample != None and refidx == -1:
            print("Error: REF sample not found", file=sys.stderr)
            sys.exit(-1)
    def ref_alt_coherence(self, refSample, altSample):
        if refSample != None and altSample == None:
            print("Error : Cannot change the REF allele without changing the ALT allele", file=sys.stderr)
            sys.exit(-1)
    def alt_integrity(self, altSample, altidx):
        if altSample != None and altidx == -1:
            print("Error : ALT sample not found", file=sys.stderr)
            sys.exit(-1)




if __name__ == "__main__":

    # Read command line arguments
    opts = get_args() 
    vcfFile = None
    refSample = None
    altSample = None
    exclusion = None
    filt_qual = 2
    verbose = False

    test_input_params = Check_input_parameters()
   
    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--vcf"):
            vcfFile = arg
        elif opt in ("-r", "--ref"):
            refSample = arg
        elif opt in ("-a", "--alt"):
            altSample = arg
        elif opt in ("-f", "--filt"):
            filt_qual = int(arg)
        elif opt in ("-x", "--exclude"):
            exclusion = arg
        elif opt in ("-v", "--verbose"):
                verbose = True
        else:
            assert False, "unhandled option"

    if verbose:        
        print("## Loading VCF file {} ...".format(vcfFile), file=sys.stderr)
        print("## Loading VCF file {} ...".format(vcfFile), file=sys.stderr)
        print("## Alt = {}".format(altSample), file=sys.stderr)
        print("## Filtering level = {}".format(filt_qual), file=sys.stderr)
    
    if filt_qual > 2:
        print("Error: --filt", file=sys.stderr)
        usage()
        sys.exit()

    if vcfFile.endswith('.gz') or vcfFile.endswith('.gzip'):
        vcf_handle = gzip.open(vcfFile)
    else:
        vcf_handle = open(vcfFile)
    header = []
    samples = []
    altidx = -1
    refidx = -1
    contaidx=[]
    
    var_counter = 0
    snp_counter = 0
    hetero_counter = 0
    badqual_counter = 0
    undefined_counter = 0
    nonspe_counter = 0
    conta_counter = 0

    for line in vcf_handle:
        try:
            line = line.rstrip()
        except TypeError:
             continue
        #print >> sys.stderr, line
        ## for now we don't care about the header 
        if line.startswith('##'):
            if refSample is not None and line.startswith("##reference="):
                print("##reference= {} [ {} ]".format(vcfFile, refSample))
            else:
                print(line)
            continue
        elif line.startswith('#'):
            header = line.split('\t')
            samples = [s.split('.')[0] for s in header[9:]]
            for  i in range(len(samples)):
                if samples[i] == refSample:
                    refidx = i
                elif samples[i] == altSample:
                    altidx = i
                elif exclusion is not None:
                    ## conta idx
                    exs = exclusion.split(",")
                    for i in range(len(exs)):
                        ct = exs[i]
                        if samples[i] == ct:
                            contaidx.append(i)
                            if verbose:
                                print("## Potential Contaminant(s) = {}".format(ct), file=sys.stderr)                           


            ## Check if Bl6 is in the conta list
            if exclusion is not None:
                exs = exclusion.split(",")

                for i in range(len(exs)):
                    ct = exs[i]
                    if ct == "REF":
                        contaidx.append(-1)
                        if verbose:
                            print("## Potential Contaminant(s) = REF", file=sys.stderr)

            ## Check input parameters
            test_input_params.ref_integrity(refSample, refidx)
            test_input_params.ref_alt_coherence(refSample, altSample)
            test_input_params.alt_integrity(altSample, altidx)
            #if refSample != None and refidx == -1:
                #print("Error : REF sample not found", file=sys.stderr)
                #sys.exit(-1)

            #if refSample != None and altSample == None:
                #print("Error : Cannot change the REF allele without changing the ALT allele", file=sys.stderr)
                #sys.exit(-1)

            #if altSample != None and altidx == -1:
                #print("Error : ALT sample not found", file=sys.stderr)
                #sys.exit(-1)

            if refidx != -1:
                print(str(' '.join(header[0:9])) + " " + refSample + "-" + altSample + "-F1")
            else:
                print(str(' '.join(header[0:9])) + " " + "REF-" + altSample + "-F1")
            continue
        else:
            if altidx == -1 :
                print("Error : ALT name not found", file=sys.stderr)
                sys.exit(-1)

            fields = line.split('\t',9)
            var_counter += 1

            ## init list of contaminant
            contg = []

            ## check chromosomes name
            if re.compile('^chr').match(fields[0]):
                chrom = fields[0]
            else:
                chrom = "chr" + str(fields[0])

            ## Filter on PASS
            if filt_qual != 2 or (filt_qual == 2 and fields[6]=="PASS"):
                ## Check format for first variant
                if var_counter == 1:
                    f = fields[8].split(':')
                    if f[0] != "GT":
                        print("Error : GT is expected to be at first index", file=sys.stderr)
                        sys.exit(-1)
                    if filt_qual == 1 and f[len(f)-1] != "FI":
                        print("Error : FI is expected to be at the last index", file=sys.stderr)
                        sys.exit(-1)
            
                genotypes  = fields[9].split('\t')
                altg = genotypes[altidx].split(':')
                altfi = altg[len(altg)-1]

                if refidx != -1:
                    refg = genotypes[refidx].split(':')
                    reffi = refg[len(refg)-1]
                else:
                    refg = ["0/0"]
                    reffi =  "1"
                    
                if len(contaidx) > 0:
                    for  i in range(len(contaidx)):
                        if contaidx[i] == -1:
                            contg.append("0/0")
                        else:
                            cg = genotypes[contaidx[i]].split(':')
                            cfi = cg[len(cg)-1]
                            if filt_qual != 1 or filt_qual == 1 and cfi == str(1):
                                contg.append(cg[0])

                ## Filter on FI field
                if filt_qual != 1 or (filt_qual == 1 and reffi == str(1) and altfi == str(1)):
                    #print "---------"
                    #print refg
                    #print altg
                    #print fields
                    geno = get_filter_snp_gt(refg[0], altg[0], fields[3], fields[4], contg)

                    if geno == -1:
                        undefined_counter += 1 
                    elif geno == -2:
                        hetero_counter += 1
                    elif geno == -3:
                        nonspe_counter += 1
                    elif geno == -4:
                        conta_counter += 1
                    else:
                        snp_counter += 1
                        #altg[0]="1/1"
                        ##print chrom + "\t" + fields[1] + "\t" + fields[2] + "\t" + geno[0] + "\t" + geno[1] + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8] + "\t" + ":".join(altg)
                        print("{}\t{}\t{}\t{}\t{}\
                              \t{}\t{}\t{}\tGT\t0|1".format(chrom, fields[1], fields[2],\
                                                            geno[0], geno[1], fields[5], fields[6], fields[7]))

                else:
                    badqual_counter += 1
            else:
                badqual_counter += 1
           
            if (verbose and var_counter % 100000 == 0):
                print("##{}".format(var_counter), file=sys.stderr)


    if verbose:
         print("## extract SNPs report", file=sys.stderr)
         print("## Total Number of SNPs ={}".format(var_counter), file=sys.stderr)
         print("## Number of reported SNPs ={}".format(snp_counter), file=sys.stderr)
         print("## -------------------------", file=sys.stderr)
         print("## Number of non discriminant SNPs ={}".format(nonspe_counter), file=sys.stderr)
         print("## Number of heterozygous SNPs ={}".format(hetero_counter), file=sys.stderr)
         print("## Number of undefined genotype SNPs ={}".format(undefined_counter), file=sys.stderr)
         print("## Number of bad quality SNPs =".format(badqual_counter), file=sys.stderr)
         print("## Number of potential contaminant SNPs ={}".format(conta_counter), file=sys.stderr)


    vcf_handle.close()
   

