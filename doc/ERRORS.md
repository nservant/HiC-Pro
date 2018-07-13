## Frequently Reported Errors

Please see the [HiC-Pro forum](https://groups.google.com/forum/#!forum/hic-pro) and the [github issue page](https://github.com/nservant/HiC-Pro).

### HiC-Pro does not generate any maps

HiC-Pro is using the chrom.sizes files to build the map.  
Be sure that your chromosome names are the same in all annotations files (bowtie2 indexes, restriction fragments, chromosome sizes, etc.)


### Error - sort: stray character in field spec: invalid field specification '2,2V'

HiC-Pro is using the option -V to sort the file per genomic location.  
On some system, the sort command does not have the -V option. In this case, please install a version of the GNU sort command which support the -V option.

### Error - pysam: AttributeError: 'module' object has no attribute 'AlignmentFile'

HiC-Pro is based on the pysam python library. Old pysam version does not support the 'AlignmentFile' class. Please, update your pysam version, with a recent one.

### Error in Paring of R1 and R2 tags

Reads from the R1 and R2 input files should be have the name.  
BAM files generated before the pairing step should be sorted by name, as the pairing is done line by line.
