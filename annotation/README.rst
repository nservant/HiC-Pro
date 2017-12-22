## How to generate restriction fragment files ?

Digest the reference genome by the provided restriction enzymes(s) and generate a BED file with the list of restriction fragments after digestion.
This file can then be used by HiC-Pro (GENOME_FRAGMENT) for the data processing.
Note that the cutting site of the restriction enzyme has to be specified using the '^' character.
The restriction enzymes HindIII, DpnII and BglII are encoded within the script and are therefore recognized if specified to the program.
Finally, multiple restriction enzymes can also be provided.

.. code-block:: bash

   ## Digest the mm9 genome by HindIII
      HICPRO_PATH/bin/utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta

   ## The same ...
      HICPRO_PATH/bin/utils/digest_genome.py -r hindiii -o mm9_hindiii.bed mm9.fasta

   ## Double digestion, HindIII + DpnII
      HICPRO_PATH/bin/utils/digest_genome.py -r hindiii dpnii -o mm9_hindiii_dpnii.bed mm9.fasta


## How to generate chromosome size files ?


The chromosome size files can be easily find on the UCSC genome browser.
Otherwise, they can also be generated using other tools such R.

require(BSgenome.Hsapiens.UCSC.hg19)

human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:25]
chrom.size <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[human_chr]
write.table(chrom_size, file="chrom_hg19.sizes", quote=FALSE, col.names=FALSE, sep="\t")
