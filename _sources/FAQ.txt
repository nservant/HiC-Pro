FAQ
===


**HiC-Pro failed. How to find what's going wrong ?**

A log file is generated when you run HiC-Pro. Its name is specified in the configuration file. In addition, some parts of the pipeline generate specific logs which are available in the *log* folder

**How can I split my .fastq files into smaller files ?**

See the `HiC-Pro Utilities <UTILS.rst>`_ which is baed on the split unix command.
For more information, see :code:`man split`

**How can I generate my annotation files ?**

HiC-Pro requires two annotation files.

* The chromosomes size are usually available through annotation website, such as the UCSC Genome Browser:

   - `hg19 <http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=13085504&chromInfoPage=>`_

   - `mm9 <http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&chromInfoPage=>`_

   - ...

Another way to generate this file is, for instance, to use the R environment.

.. code-block:: guess

   ###
   ## How to generate chromosome size files ?
   ### 
   require(BSgenome.Hsapiens.UCSC.hg19)

   human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
   chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[human_chr]
   write.table(chrom_size, file="chrom_hg19.sizes", quote=FALSE, col.names=FALSE, sep="\t")


* The restriction fragments file has to be generated according to the reference genome, and the restriction enzyme(s) used to generate the Hi-C data.

Here is the way to generate such file using the `HiTC <http://bioconductor.org/packages/release/bioc/html/HiTC.html>`_ BioConductor package.

.. code-block:: guess

   ###
   ## How to generate restriction fragment files ?
   ### 

   require(HiTC)
   require(rtracklayer)
   require(BSgenome.Hsapiens.UCSC.hg19)

   human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
   resFrag <- getRestrictionFragmentsPerChromosome(resSite="AAGCTT", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg19")
   allRF <- do.call("c",resFrag)
   names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
   export(allRF, format="bed", con="HindIII_resfrag_hg19.bed")


**Why HiC-Pro need to be run in two steps in parallel mode ?**

Th HiC-Pro pipeline is divided into two main steps. The first part of the pipeline is 'fastq' based, meaning that the same anlaysis will be performed for all fastq files.
This part can be easily parallelized per fastq, with at the end, a list of valid interactions per fastq file.
The second step of the pipeline is 'sample' based. All lists of valid interactions from the same sample are merged in order to build and normalize the maps.
At that time, this second step is not time consuming, and we do not parallelize it although a per sample parallelization migth be a good idea.
So, because these two steps as either 'fastq' based or 'sample' based, we need to separate them during the parallele processing.

**Does HiC-Pro support other mapper ?**

No. HiC-Pro is only based on the bowtie2 mapper.
However, note that HiC-Pro can be run from aligned data. In this case, the input path (-i) must be a BAM folder, and the analysis has to be run step-by-step.

**How can I create N-masked genome for allele-specific analysis ?**

The allele specific mode of HiC-Pro is based on a N-masked genome. Meaning that all SNPs information which can be use to distinguish parental haplotypes have to be masked. This masking can be performed in 3 steps:
1. Extract relevant SNPs information. See the `extract_snps.py <doc/UTILS.rst>`_ utility for Mouse Sanger data. For Human data, you can use phasing data, or SNPs information available from public ressources, as the `Illumina Platinum Project <http://www.illumina.com/platinumgenomes/>`_, the `1K Genome Project <http://www.1000genomes.org/>`_ or the `GATK resource bundle <https://www.broadinstitute.org/gatk/guide/article.php?id=1215>`_.
2. Mask the fasta genome. To do so, simply use the bedtools `maskfasta <http://bedtools.readthedocs.org/en/latest/content/tools/maskfasta.html>`_ utility.
3. Then, create your bowtie2 indexes from the masked fasta file.
