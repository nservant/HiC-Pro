## HiC-Pro Utilities

HiC-Pro provides a couple of utilities which are not part of the main pipeline, but was designed to help the user to perform a couple of tasks.  
All utilities are installed in /bin/utils/.


### split_reads.py OR 'How can I split my reads in chunck ?'

Split reads is highly recommanded in HiC-pro parallel mode.  
To split the reads in chunck, you can simply used the split_reads.py utility.

```
   ## split the fastq file in blocks of 10M reads
   HICPRO_PATH/bin/utils/split_reads.py --results_folder OUTPUT --nreads READS_NB INPUT_FASTQ
```

### extract_snps.py OR 'How can I extract SNPs information from phasing data ?'

This script was designed to extract SNPs information from a VCF file.  
As a example, the Sanger database provides a list of SNPs which can be used to differentiate several Mouse strains (see ftp://ftp-mouse.sanger.ac.uk/current_snps/).  
The single VCF contains all SNPs information for all Mouse strains. The default reference allele is the C57 black strain.  
This script is able to extract informative and high quality SNPs (-f option) which can be then used in allele specific analysis to distinguish parental origin.

```
   ## Extract SNPs information for CASTEiJ/129S1 cross
   HICPRO_PATH/bin/utils/extract_snps.py -i mgp.v2.snps.annot.reformat.vcf -r CASTEij -a 129S1 > snps_CASTEiJ_129S1.vcf

   ## Extract SNPs information for C57_b6/FVB_NJ
   HICPRO_PATH/bin/utils/extract_snps.py -i mgp.v2.snps.annot.reformat.vcf -a FVB_NJ > snps_C57b6_FVBNJ.vcf
```

### digest_genome.py OR 'How can I generate the list of restriction fragments after genome digestion ?'

Digest the reference genome by the provided restriction enzymes(s) and generate a BED file with the list of restriction fragments after digestion.  
This file can then be used by HiC-Pro (GENOME_FRAGMENT) for the data processing.  
Note that the cutting site of the restriction enzyme has to be specified using the '^' character.  
The restriction enzymes HindIII, DpnII, BglII and MboI are encoded within the script and are therefore recognized if specified to the program.  
Finally, note that multiple restriction enzymes can also be provided (space separated), and that since v2.11.3, 'N' based are recognized and replaced by A,C,G,T  

```
   ## Digest the mm9 genome by HindIII
   HICPRO_PATH/bin/utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta

   ## The same ...
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii -o mm9_hindiii.bed mm9.fasta

   ## Double digestion, HindIII + DpnII
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii dpnii -o mm9_hindiii_dpnii.bed mm9.fasta
```


### make_viewpoints.py OR 'How can I generate a BED profile from a given viewpoints ?'

The make_viewpoints.py script was initially designed in the context of capture-C data.  
For a list of anchors (i.e. capture sites), it allows to generate track files with all interacting fragments.  
Regions around the capture are usually excluded from the profile ('-e' parameter)

```
   ## Generate a viewpoints from capture site
   HICPRO_PATH/bin/utils/make_viewpoints -i hicpro_res/hic_results/data/dixon_2M/dixon_2M_allValidPairs -f HICPRO_PATH/data_info/HindIII_resfrag_hg19.bed -t mycapture.bed -e 1000 -o OUTPUT_DIR/output_ -v
```

### hicpro2juicebox.sh OR 'How can I load my HiC-Pro data into Juicebox visualization software ?'

The hicpro2juicebox.sh utility allows to generate input file for Juicebox.  
It can be use both for restriction fragment Hi-C or Dnase Hi-C by specifying the -f FRAGMENT_FILE option. Note that in this case, it is advice to use the HiC-Pro annotation file, as the fragment name is expected to be HiC_CHROMOSOME_FRAGMENTNUMBER.  
This utility requires HiC-Pro version 2.7.6 or later, and the installation of [Juicer command line tools](https://github.com/theaidenlab/juicer/wiki/Download)

```
   ## Convert HiC-Pro output to Juicebox input
   HICPRO_PATH/bin/utils/hicpro2juicebox.sh -i hicpro_res/hic_results/data/dixon_2M/dixon_2M_allValidPairs -g chrom_hg19.sizes -j /usr/local/juicer/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar

   ## Convert HiC-Pro output to Juicebox input up to restriction fragment resolution
   HICPRO_PATH/bin/utils/hicpro2juicebox.sh -i hicpro_res/hic_results/data/dixon_2M/dixon_2M_allValidPairs -g chrom_hg19.sizes -j /usr/local/juicer/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar -f  HICPRO_PATH/data_info/HindIII_resfrag_hg19.bed
```

### sparseToDense.py OR 'How can I convert HiC-Pro output into dense format ?'

ALlows to convert data in sparse symmetric format into dense matrices. This convertion can be useful for downstream analysis such as TADs calling using the directionaly index method (Dixon et al. 2012). The utility can also be used to extract intra-chromosomal maps at dense format.

```
  ## Convert to dense format
  HICPRO_PATH/bin/utils/sparseToDense.py -b hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed hic_results/matrix/dixon_2M/iced/1000000/dixon_2M_1000000_iced.matrix

  ## Convert todense format per chromosome
  HICPRO_PATH/bin/utils/sparseToDense.py -b hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed hic_results/matrix/dixon_2M/iced/1000000/dixon_2M_1000000_iced.matrix --perchr

  ## Convert into TADs caller input from Dixon et al.
  HICPRO_PATH/bin/utils/sparseToDense.py -b hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed hic_results/matrix/dixon_2M/iced/1000000/dixon_2M_1000000_iced.matrix --perchr --di
```

### hicpro2fithic.py OR 'How can I use Fit-Hi-C after HiC-Pro processing ?'

Convert HiC-Pro output to Fit-Hi-C input (Ay et al. 2014)

```
  ## Whith IC bias vector
  HICPRO_PATH/bin/utils/hicpro2fithic.py -i hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000.matrix -b hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000_abs.bed -s hic_results/matrix/dixon_2M/iced/1000000/dixon_2M_1000000_iced.matrix.biases
```

### hicpro2higlass.sh OR 'How can I load my HiC-Pro output into the Higlass visualization tool ?'

First be sure that Higlass is install on the environment and that the cooler python package is available.  
See http://gehlenborglab.org/research/projects/higlass/ for higlass installation.  
See https://github.com/mirnylab/cooler for .cool Hi-C data format. The path to cooler utility must be defined in your PATH.

```
   ## Convert matrix file into .cool file
   HICPRO_PATH/bin/utils/hicpro2higlass.sh -i hic_results/matrix/dixon_2M/raw/1000000/dixon_2M_1000000.matrix -r 1000000 -c chrom_hg19.sizes -n

   ## Convert matrix file into .cool file
   HICPRO_PATH/bin/utils/hicpro2higlass.sh -i hic_results/data/dixon_2M/dixon_2M_allValidPairs -r 40000 -c chrom_hg19.sizes -n
```

Note that you should use the `--float` option for iced normalized data.

### split_sparse.py OR 'How can I split my HiC-Pro sparse format ?'

This utility allows to split a genome-wide contact maps generated by HiC-Pro into per-chromosome interaction maps.
This can be especially useful for TADs calling.

```
HICPRO_PATH/bin/utils/split_sparse.py -b mydata_abs.bed mydata.matrix
```
