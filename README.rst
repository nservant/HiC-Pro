.. _QS:

.. Nicolas Servant
.. HiC-Pro
.. v2.5.0
.. 15-04-02

HiC-Pro Quick Start Guide
*************************

This page is just a quick start guide, please read the full `online manual <http://nservant.github.io/HiC-Pro/>`_ for more information !

See NEWS for information about changes in this and previous versions

See LOGBOOK for details about the HiC-Pro developments

See LICENSE for license information


What is HiC-Pro ?
=================

| HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. Since version 2.7.0, HiC-Pro can analyse data from digestion protocols as well as data from protocols that do not require restriction enzyme such as DNase Hi-C.
| The pipeline is flexible, scalable and optimized. It can operate either on a single laptop or on a computational cluster using the PBS-Torque scheduler. HiC-Pro is sequential and each step of the workflow can be run independantly.
| HiC-Pro includes a fast implementatation of the ICE normalization method (see the `iced <https://github.com/hiclib/iced>`_ python library for more information).
| In addition, HiC-Pro can use phasing data to build :ref:`allele specific contact maps <AS>`.

If you use HiC-Pro, please cite :

HiC-Pro: An optimized and flexible pipeline for Hi-C processing. *Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E. Genome Biology 2015, 16:259* doi:10.1186/s13059-015-0831-x
`http://www.genomebiology.com/2015/16/1/259 <http://www.genomebiology.com/2015/16/1/259>`_

For any question about HiC-Pro, please contact nicolas.servant@curie.fr

How to install it ?
===================

The HiC-Pro pipeline requires the following dependencies :

* The `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ mapper
* Python (>2.7) with *pysam*, *bx*, *numpy*, and *scipy* libraries
* R with the *RColorBrewer* and *ggplot2* packages
* g++ compiler
* Samtools (>0.1.19)

Bowtie >2.2.2 is strongly recommanded for allele specific analysis.
| To install HiC-Pro:

.. code-block:: guess

  tar -zxvf HiC-Pro-master.tar.gz
  cd HiC-Pro-master
  make CONFIG_SYS=config-install.txt install

| Note that if some of these dependencies are not installed (i.e. not detected in the $PATH), HiC-Pro will try to install them.
| You can also edit the *config-install.txt* file and manually defined the paths to dependencies.

+---------------+------------------------------------------------------------+
| SYSTEM CONFIGURATION                                                       |
+===============+============================================================+
| PREFIX        | Installation path                                          |
+---------------+------------------------------------------------------------+
| BOWTIE2_PATH  | Full path the bowtie2 installation directory               |
+---------------+------------------------------------------------------------+
| SAMTOOLS_PATH | Full path to the samtools installation directory (>0.1.19) |
+---------------+------------------------------------------------------------+
| R_PATH        | Full path to the R installation directory                  |
+---------------+------------------------------------------------------------+
| PYTHON_PATH   | Full path to the python installation directory (>2.7)      |
+---------------+------------------------------------------------------------+


Annotation Files
================

In order to process the raw data, HiC-Pro requires three annotation files. Note that the pipeline is provided with some Human and Mouse annotation files.

1. **A BED file** of the restriction fragments after digestion. This file depends both of the restriction enzyme and the reference genome. See the :ref:`FAQ <FAQ>` and the :ref:`HiC-Pro utilities <UTILS>` for details about how to generate this file. A few annotation files are provided with the HiC-Pro sources.

::

   chr1   0       16007   HIC_chr1_1    0   +
   chr1   16007   24571   HIC_chr1_2    0   +
   chr1   24571   27981   HIC_chr1_3    0   +
   chr1   27981   30429   HIC_chr1_4    0   +
   chr1   30429   32153   HIC_chr1_5    0   +
   chr1   32153   32774   HIC_chr1_6    0   +
   chr1   32774   37752   HIC_chr1_7    0   +
   chr1   37752   38369   HIC_chr1_8    0   +
   chr1   38369   38791   HIC_chr1_9    0   +
   chr1   38791   39255   HIC_chr1_10   0   +
   (...)

2. **A table file** of chromosomes' size.

::

   chr1    249250621
   chr2    243199373
   chr3    198022430
   chr4    191154276
   chr5    180915260
   chr6    171115067
   chr7    159138663
   chr8    146364022
   chr9    141213431
   chr10   135534747
   (...)

3. **The bowtie2 indexes**. See `the bowtie2 manual page <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ for details about how to create such indexes.

How to use it ?
===============

0. First have a look at the help message !

.. code-block:: guess

  HiC-Pro --help
  usage : HiC-Pro -i INPUT -o OUTPUT -c CONFIG [-s ANALYSIS_STEP] [-p] [-h] [-v]
  Use option -h|--help for more information

  HiC-Pro 2.7.0
  ---------------
  OPTIONS

   -i|--input INPUT : input data folder; Must contains a folder per sample with input files
   -o|--output OUTPUT : output folder
   -c|--conf CONFIG : configuration file for Hi-C processing
   [-p|--parallel] : if specified run HiC-Pro in PBS/Torque mode
   [-s|--step ANALYSIS_STEP] : run only a subset of the HiC-Pro workflow; if not specified the complete workflow is run
      mapping: perform reads alignment
      proc_hic: perform Hi-C filtering
      quality_checks: run Hi-C quality control plots
      build_contact_maps: build raw inter/intrachromosomal contact maps
      ice_norm: run ICE normalization on contact maps
   [-h|--help]: help
   [-v|--version]: version

1. Copy and edit the configuration file *'config-hicpro.txt'* in your local folder. See the :ref:`manual <MANUAL>` for details about the configuration file
2. Put all input files in a rawdata folder. The input files have to be organized with a folder per sample.
3. Run HiC-Pro

* **On your laptop**

.. code-block:: guess

    MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE


* **Using a cluster (PBS)**

.. code-block:: guess

   MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -p



You will get the following message :

.. code-block:: guess

  Please run HiC-Pro in two steps :
  1- The following command will launch the parallel workflow through 12 torque jobs:
  qsub HiCPro_step1.sh
  2- The second command will merge all outputs to generate the contact maps:
  qsub HiCPro_step2.sh


Execute the displayed command from the output folder:

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
   normContactMaps.sh -c /bioinfo/users/nservant/projects_dev/HiC-Pro/config_test.txt >> hicpro_IRM90_rep1_split.log 

   real	3m23.902s
   user	5m22.956s
   sys	0m40.243s

   ## All results are available in HiC_Pro_v2.4.0_test

