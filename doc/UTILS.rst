..  _UTILS:

HiC-Pro Utilities
=================


| HiC-Pro provides a couple of utilities which are not part of the main pipeline, but was designed to help the user to perform a couple of tasks.
| All utilities are installed in /bin/utils/.


1- split_reads.py
-----------------

or **How can I split my reads in chunck ?**

| Split reads is highly recommanded in HiC-pro parallel mode.
| To split the reads in chunck, you can simply used the split_reads.py utility.

.. code-block:: bash

   ## split the fastq file in blocks of 10M reads
   HICPRO_PATH/bin/utils/split_reads.py --results_folder OUTPUT --nreads READS_NB INPUT_FASTQ


2- extract_snps.py
------------------
or **How can I extract SNPs information from phasing data ?**

| This script was designed to extract SNPs information from a VCF file.
| As a example, the Sanger database provides a list of SNPs which can be used to differentiate several Mouse strains (see ftp://ftp-mouse.sanger.ac.uk/current_snps/).
| The single VCF contains all SNPs information for all Mouse strains. The default reference allele is the C57 black strain.
| This script is able to extract informative and high quality SNPs (-f option) which can be then used in allele specific analysis to distinguish parental origin.

.. code-block:: bash

   ## Extract SNPs information for CASTEiJ/129S1 cross
   HICPRO_PATH/bin/utils/split_reads.py -i mgp.v2.snps.annot.reformat.vcf -r CASTEij -a 129S1 > snps_CASTEiJ_129S1.vcf

   ## Extract SNPs information for C57_b6/FVB_NJ
   HICPRO_PATH/bin/utils/split_reads.py -i mgp.v2.snps.annot.reformat.vcf -a FVB_NJ > snps_CASTEiJ_129S1.vcf



2- digest_genome.py
-------------------
or **How can I generate the list of restriction fragments after genome digestion ?**

| Digest the reference genome by the provided restriction enzymes(s) and generate a BED file with the list of restriction fragments after digestion.
| This file can then be used by HiC-Pro (GENOME_FRAGMENT) for the data processing.
| Note that the cutting site of the restriction enzyme has to be specified using the '^' character.
| The restriction enzymes HindIII, DpnII and BglII are encoded within the script and are therefore recognized if specified to the program.
| Finally, multiple restriction enzymes can also be provided.

.. code-block:: bash

   ## Digest the mm9 genome by HindIII
   HICPRO_PATH/bin/utils/digest_genome.py -r A^AGCTT -o mm9_hindiii.bed mm9.fasta

   ## The same ...
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii -o mm9_hindiii.bed mm9.fasta

   ## Double digestion, HindIII + DpnII
   HICPRO_PATH/bin/utils/digest_genome.py -r hindiii dpnii -o mm9_hindiii_dpnii.bed mm9.fasta




