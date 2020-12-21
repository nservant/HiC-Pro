#! /usr/bin/env python

import argparse
import math
import os
import gzip

# Created by Arya Kaul - 1/12/2017
# Modified by Ferhat Ay - 1/13/2017
# Modified by Nicolas Servant - 1/23/2017
# Modified by Ferhat Ay - 6/5/2017 - added resolution (-r <res>) argument to avoid some problems with inferring it from the first entry of bedFile
# Modified by Ferhat Ay - 6/7/2017 - hitCount to ints, and gzip output 
# Modified by Arya Kaul - 12/21/2017 - bug fix when no -o option specified 


def outputfithicform(bedPath, matrixPath, intCPath, fragMapPath, biasVectorPath=None, biasVectorOutput=None, res=0):
	print("Loading matrix file...")
	fragDic = {}
	# resolution of data to be determined if res=0 at this point
	with open(bedPath, 'r') as bedFile:
		for lines in bedFile:
			line = lines.rstrip().split()
			chrNum, start, en = line[0], line[1], line[2]
			if res == 0: res = int(en) - int(start) # edge case but if first fragment in file is smaller than res, there is a problem
			mid = int(start) + int(res/2)
			index = int(line[3])
			fragDic[index] = [chrNum, start, mid, 0] # last field is total contact count: tcc

	lineCount = 0
	with open(matrixPath, 'r') as matrixFile:
		with gzip.open(intCPath, 'wb') as interactionCountsFile:
			for lines in matrixFile:
				line = lines.rstrip().split()
				i, j = int(line[0]), int(line[1])
				cc = float(line[2]) # this can be float or int
				fragDic[i][3] += cc
				fragDic[j][3] += cc
				if cc == int(cc):
					cc = int(cc) # make sure to convert to integer if it ends with ".0"
				interactionCountsFile.write("{}\t{}\t{}\t{}\t{}\n".format(str(fragDic[i][0]),
                                                                          str(fragDic[i][2]),
                                                                          str(fragDic[j][0]),
                                                                          str(fragDic[j][2]),
                                                                          str(cc)).encode())
				lineCount += 1
				if lineCount % 1000000 == 0:
					print("{} million lines read".format(int(lineCount/1000000)))

	with gzip.open(fragMapPath, 'wb') as fragmentMappabilityFile:
		for indices in sorted(fragDic): # sorted so that we start with the smallest index
			toWrite = 0
			if fragDic[indices][3] > 0:
				toWrite = 1
			fragmentMappabilityFile.write("{}\t{}\t{}\t{}\t{}\n".format(fragDic[indices][0],
                                                                        fragDic[indices][1],
                                                                        fragDic[indices][2],
                                                                        int(fragDic[indices][3]),
                                                                        toWrite).encode()) 

	if biasVectorPath is not None and biasVectorOutput is not None:
		print("Converting bias file...")

		biasVec = [0,0] # bias sum and biasCount
		biasDic = {} # bias for each index
		i = 1 # one-based indices
		with open(biasVectorPath, 'r') as biasVectorFile:
				for lines in biasVectorFile:
					value = float(lines.rstrip()) #just one entry that can be nan or a float
					index = int(i)
					i += 1
					biasDic[index] = value
					if not math.isnan(value):
						biasVec[0] += value #sum
						biasVec[1] += 1 # count
	#

		# Centering the bias values on 1.
		biasAvg=biasVec[0] / biasVec[1]

		with gzip.open(biasVectorOutput, 'wb') as biasVectorOutputFile:
			for index in sorted(biasDic):
				value = biasDic[index]
				if not math.isnan(value):
					value = value/biasAvg
				else: 
					value = -1
				biasVectorOutputFile.write("{}\t{}\t{}\n".format(str(fragDic[index][0]),
                                                                 str(fragDic[index][2]),
                                                                 str(value)+'\n').encode())
	print("Conversion from HiC-Pro to Fit-Hi-C format completed")

#outputfithicform(args.bedPath, args.matrixPath, args.intCPath, args.fragMapPath, args.biasVectorPathandOutput[0], args.biasVectorPathandOutput[1])

def main():
	# Example without bias files
	outputfithicform('raw/1000000/hIMR90_HindIII_r1_1000000_abs.bed', 'raw/1000000/hIMR90_HindIII_r1_1000000.matrix', 'fithic.interactionCounts', 'fithic.fragmentMappability')
	# Example with bias files
	outputfithicform('raw/1000000/hIMR90_HindIII_r1_1000000_abs.bed', 'raw/1000000/hIMR90_HindIII_r1_1000000.matrix', 'fithic.interactionCounts', 'fithic.fragmentMappability','hicpro.biases','fithic.biases')

if __name__=="__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--matrix", help="Input matrix file with raw contact frequencies.", required=True)
	parser.add_argument("-b", "--bed", help="BED file with bins coordinates.", required=True)
	parser.add_argument("-s", "--bias", help="The bias file provided after IC normalization.", default=None)
	parser.add_argument("-o", "--output", help="Output path", default=".")
	parser.add_argument("-r", "--resolution", help="Resolution of the matrix", type=int, default=0) # 0 means it is inferred from fragments file's first entry

	args = parser.parse_args()

	icounts_output = os.path.join(args.output + "/fithic.interactionCounts.gz")
	fragmap_output = os.path.join(args.output + "/fithic.fragmentMappability.gz")
	bias_output = None

	if args.bias is not None:
		bias_output = os.path.join(args.output + "/fithic.biases.gz")

	outputfithicform(args.bed, args.matrix, icounts_output, fragmap_output, args.bias, bias_output, args.resolution)

