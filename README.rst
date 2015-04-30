.. Nicolas Servant
.. HiC-Pro
.. v2.5.0
.. 15-04-02

HiC-Pro Quick Start Guide
*************************

See NEWS for information about changes in this and previous versions

See LOGBOOK for details about the HiC-Pro developmens

See LICENSE for license information

What is HiC-Pro ?
=================

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. 
The pipeline is flexible, scalable and optimized. It can operate either on a single laptop or on a computational cluster using the PBS-Torque scheduler

If you use HiC-Pro, please cite :

HiC-Pro: An optimized and flexible pipeline for Hi-C processing. *Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E.*. 2015. submitted

How to install it ?
===================

The HiC-Pro pipeline requires the following dependencies :

* The bowtie2 mapper (or any other mapper)
* Python with *pysam*, *bx*, *numpy*, and *scipy* libraries
* R with the *RColorBrewer* and *ggplot2* packages
* g++ compiler
* Samtools (>0.1.18)

To install HiC-Pro:

.. code-block:: guess

  tar -zxvf HiC-Pro.2.3.1.tar.gz
  cd HiC-Pro_2.3.1
  make CONFIG_SYS=config-install.txt install



Note that if some of these dependencies are not installed (i.e. not detected in the $PATH), HiC-Pro will try to install them.
You can also edit the *config-install.txt* file and manually defined the paths to dependencies.

Annotation Files
================

In order to process the raw data, HiC-Pro requires three annotation files :

1. A BED file of the restriction fragments after digestion. This file depends both of the restriction enzyme and the reference genome. See the `FAQ <../html/FAQ.html>`_ for details about how to generate this file. A few annotation files are provided with the HiC-Pro sources.
2. A table file of chromosomes' size.
3. The bowtie2 indexes. See `the bowtie2 manual page <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ for details about how to create such indexes.



How to use it ?
===============

1. Copy and edit the configuration file *'config-hicpro.txt'* in your local folder.
2. Put all fastq files in a rawdata folder. Each fastq file has to be put in a folder per sample.
3. Run HiC-Pro

  * Without PBS-Torque

  .. code-block:: guess

    MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE
  
  * With PBS-Torque

  .. code-block:: guess

   MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -p



You will get the following message :

.. code-block:: guess

  Please run HiC-Pro in two steps :
  1- The following command will launch the parallel workflow through 12 torque jobs:
  qsub HiCPro_step1.sh
  2- The second command will merge all outputs to generate the contact maps:
  qsub HiCPro_step2.sh


Execute the displayed command:

.. code-block:: guess

  qsub HiCPro_step1.sh
  774410[].torque.curie.fr


Then wait for the torque mails... :)
Once executed succesfully (may take several hours), then type:

.. code-block:: guess

  qsub HiCPro_step2.sh


Test Dataset
============

Small fastq files (2M reads) extracted from the Dixon et al. 2012 paper are available for test.

.. code-block:: guess

   ## Get the data. Will download a test_data folder and a configuration file
   wget -r -nH --cut-dirs=2 --no-parent --reject="index.html*" https://zerkalo.curie.fr/partage/HiC-Pro/

   ## Run HiC-Pro
   time HICPRO_INSTALL_DIR/bin/HiC-Pro -i test_data -o HiC_Pro_v2.4.0_test -c config_test.txt
   Run HiC-Pro
   --------------------------------------------
   lundi 2 mars 2015, 17:00:36 (UTC+0100)
   Bowtie2 global alignment ...
   bowtie_wrap.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt -u >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:01:25 (UTC+0100)
   Bowtie2 local alignment ...
   bowtie_wrap.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt -l >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:01:41 (UTC+0100)
   Combine both alignment ...
   bowtie_combine.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:01:52 (UTC+0100)
   Bowtie2 mapping statistics for R1 and R2 tags ...
   mapping_stat.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:01:53 (UTC+0100)
   Pairing of R1 and R2 tags ...
   bowtie_pairing.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:02:22 (UTC+0100)
   Assign alignments to HindIII sites ...
   mapped_2hic_fragments.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:03:49 (UTC+0100)
   Merge multiple files from the same sample ...
   merge_valid_interactions.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:03:49 (UTC+0100)
   Make plots per sample ...
   make_plots.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log
   --------------------------------------------
   lundi 2 mars 2015, 17:03:55 (UTC+0100)
   Generate binned matrix files ...
   build_raw_maps.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt 2> logs/build_raw_maps.log
   --------------------------------------------
   lundi 2 mars 2015, 17:03:57 (UTC+0100)
   Run ICE Normalization ...
   normContactMaps.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log # 2> logs/normICE.log

   real	3m23.902s
   user	5m22.956s
   sys	0m40.243s

   ## All results are available in HiC_Pro_v2.4.0_test

