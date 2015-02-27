FAQ
===


**1- HiC-Pro failed. How to find what's going wrong ?**

A log file is generated when you run HiC-Pro. Its name is specified in the configuration file. In addition, some parts of the pipeline generate specific logs which are available in the *log* folder

**2- How can I split my .fastq files into smaller files ?**

In linux system, you can simply use the :code:`split` command.


.. code-block:: bash

   ## split the fastq file in blocks of 10M reads
   split -l 10000000 -d FILE SUFFIX


For more information, see :code:`man split`

**3- How can I generate my annotation files ?**

HiC-Pro requires two annotation files.

* The chromosomes size are usually available through annotation website, such as the UCSC Genome Browser:

   - `hg19<http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgsid=13085504&chromInfoPage=>`_

   - `mm9<http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&chromInfoPage=>`_

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
Here is the way to generate such file using the `HiTC<http://bioconductor.org/packages/release/bioc/html/HiTC.html>`_ BioConductor package.

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


