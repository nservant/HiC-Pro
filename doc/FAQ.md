## FAQ

### Can I use HiC-Pro with the Arima Hi-C Kit ?

Yes, the LIGATION_SITE variable from the configuration file can contain multiple ligation motifs (comma separated).

Regarding the Arima kit, as there is two restiction motifs (^GATC and G^ANTC), the possible ligation motifs should be GATCGATC,GANTGATC,GANTANTC,GATCANTC.

Leading to 25 possible ligation motifs when the 'N' are replaced by 'A', 'C', 'G', 'T'.
Since HiC-Pro version 2.11.3, 'N' base are automatically recognized.
For older version, please specify all possible ligation motifs.

### HiC-Pro failed. How to find what's going wrong ?

A log file is generated when you run HiC-Pro. Its name is specified in the configuration file.   
In addition, some parts of the pipeline generate specific logs which are available in the *log* folder.  
Since v2.11.0, logs have be written to be as clear as possible. Each HiC-Pro step has its own log file.


### How can I split my .fastq files into smaller files ?

See the [HiC-Pro split_reads.py utility](UTILS.md) which is based on the *split* unix command.

### How can I used my own annotations files

Simply put the full path of your annotations in the configuration file. By default HiC-Pro will check if the file exists, and if not, will look for them in its own annotation folder.

### How can I generate my annotation files ?

HiC-Pro requires two annotation files.

* The chromosomes size are usually available through annotation website, such as the UCSC Genome Browser:

   - [hg19](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=13085504&chromInfoPage=)
   - [mm9](http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&chromInfoPage=)
   - ...

Another way to generate this file is, for instance, to use the R environment.

```
   ###
   ## How to generate chromosome size files ?
   ###
   require(BSgenome.Hsapiens.UCSC.hg19)

   human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
   chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[human_chr]
   write.table(chrom_size, file="chrom_hg19.sizes", quote=FALSE, col.names=FALSE, sep="\t")
```

* The restriction fragments file has to be generated according to the reference genome, and the restriction enzyme(s) used to generate the Hi-C data.

Since version 2.7.0, HiC-Pro proposes the utility *digest_genome.py* to generate this file using as input, the fasta file and the name(s) or sequence(s) of the restriction enzyme(s).  
See the [HiC-Pro Utilities](UTILS.md) section for more details.  

Another way is to generate the list of restriction fragments is to use the [HiTC](http://bioconductor.org/packages/release/bioc/html/HiTC.html) BioConductor package.
Note that this method only works for the genomes which are already available in BioConductor and for one restriction enzyme.

The packages [HiTC](http://bioconductor.org/packages/release/bioc/html/HiTC.html), [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and [BSgenome.Hsapiens.UCSC.hg19](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg19.html) have to be installed to run the following example.

```
   ###
   ## How to generate restriction fragment files with HiTC ?
   ###

   require(HiTC)
   require(rtracklayer)
   require(BSgenome.Hsapiens.UCSC.hg19)

   human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
   resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg19")
   allRF <- do.call("c",resFrag)
   names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
   export(allRF, format="bed", con="HindIII_resfrag_hg19.bed")
```

### Why does HiC-Pro need to be run in two steps in parallel mode ?

Th HiC-Pro pipeline is divided into two main steps. The first part of the pipeline is 'fastq' based, meaning that the same anlaysis will be performed for all fastq files.  
This part can be easily parallelized per fastq, with at the end, a list of valid interactions per fastq file.  
The second step of the pipeline is 'sample' based. All lists of valid interactions from the same sample are merged in order to build and normalize the maps.  
At that time, this second step is not time consuming, and we do not parallelize it although a per sample parallelization migth be a good idea.  
So, because these two steps as either 'fastq' based or 'sample' based, we need to separate them during the parallele processing.  


### Does HiC-Pro support other mappers ?

No. HiC-Pro is only based on the bowtie2 mapper.  
However, note that HiC-Pro can be run from aligned data. In this case, the input path (-i) must be a BAM folder, and the analysis has to be run step-by-step.


### How can I create N-masked genome for allele-specific analysis ?

The allele specific mode of HiC-Pro is based on a N-masked genome. Meaning that all SNPs information which can be use to distinguish parental haplotypes have to be masked. This masking can be performed in 3 steps:

1. Extract relevant SNPs information. See the [extract_snps.py ](UTILS.md) utility for Mouse Sanger data. For Human data, you can use phasing data, or SNPs information available from public ressources, as the [Illumina Platinum Project](http://www.illumina.com/platinumgenomes/), the [1K Genome Project](http://www.1000genomes.org/) or the [GATK resource bundle](https://www.broadinstitute.org/gatk/guide/article.php?id=1215).

2. Mask the fasta genome. To do so, simply use the bedtools [maskfasta](http://bedtools.readthedocs.org/en/latest/content/tools/maskfasta.html) utility.

3. Then, create your bowtie2 indexes from the masked fasta file.


### What can I do once I have the iced contact maps

The matrix format is a standard sparse triplet format which can easily be loaded in R or matlab environment.  
For instance, the matrix can be easily loaded in the R environment using the [HiTC Bioconductor package](http://bioconductor.org/packages/release/bioc/html/HiTC.html).

```
   require(HiTC)
   ## Load Hi-C data
   x<-importC("mydata.matrix", xgi.bed="mydata_abs.bed")
   show(x)
   ## Plot X intra-chromosomal map
   mapC(HTClist(x$chrXchrX), trim.range=.9)
```

### How can I visualize the maps generated by HiC-Pro

The maps can be plotted into througth R environment. See the [HiTC Bioconductor package](http://bioconductor.org/packages/release/bioc/html/HiTC.html) and the previous question.  
The HiC-pro results are also compatible with the HiCPlotter software ([Akdemir et al. 2015](http://www.genomebiology.com/2015/16/1/198)).  
The source of HiCPlotter are available on [github](https://github.com/kcakdemir/HiCPlotter).  
Here is a small example of how to use HiCPlotter.

```
   ## Plot the genome-wide map at 1Mb resolution
   python HiCPlotter.py -f hic_results/matrix/sample1/iced/1000000/sample1_1000000_iced.matrix -o Examplegw -r 1000000 -tri 1 -bed hic_results/matrix/sample1/raw/1000000/sample1_1000000_ord.bed -n hES -wg 1 -chr chrX

   ## Plot the chrX at 150Kb resolution
   python HiCPlotter.py -f hic_results/matrix/sample1/iced/150000/sample1_150000_iced.matrix -o Exemple -r 150000 -tri 1 -bed hic_results/matrix/sample1/raw/150000/sample1_150000_ord.bed -n Test -chr chrX -ptr 1
```

Since version 2.7.6 HiC-Pro is compatible with the Juicebox viewer. See the [hicpro2juicebox utility](UTILS.md) to generate Juicebox input file from the list of valid interactions.
Since version 2.10.0 HiC-Pro is also compatible with the Higlass viewer. See the [hicpro2higlass utility](UTILS.md) to generate .cool input file.


### How much disk space is require for running HiC-Pro

In average, HiC-Pro requires 4 times more space than the raw data (more if the fastq files are compressed).  
The main part of this disk space is required for the mapping only, as HiC-Pro will map R1 and R2 independantly, in two steps. And then, merge the R1 and R2 alignments at the pairing step. So for 1 sample (R1 + R2 input files), 8 intermediate BAM files are generated.  
The final matrix files are not that big, and the triplet format used by HiC-Pro is very efficient.  
Therefore, a good habit would be to remove mapping files after running HiC-Pro.  
