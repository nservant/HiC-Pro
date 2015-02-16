.. Nicolas Servant
.. HiC-Pro
.. v2.3.1
.. 15-15-02

HiC-Pro Manual
******************

Hi-C Overview
=============


What is HiC-Pro ?
=================

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. 
The pipeline is flexible, scalable and optimized. It can operate either on a single laptop or on a computational cluster using the PBS-Torque scheduler

ADD HiC-Pro Wklfow

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
  make install



Note that if some of these dependencies are not installed (i.e. not detected in the $PATH), HiC-Pro will try to install them.
You can also edit the *config-install.txt* file and manually defined the paths to dependencies.

Annotation Files
================


How to use it ?
===============

1. Copy and edit the configuration file *'config-hicpro.txt'* in your local folder

TODO : details of the configuration file

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


How does HiC-Pro work ?
=======================

1. Reads Mapping
2. Fragment assignment and filtering
3. Map builder
4. ICE normalization

Data format
===========










