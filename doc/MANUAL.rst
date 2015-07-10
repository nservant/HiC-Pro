.. Nicolas Servant
.. HiC-Pro
.. v2.3.1
.. 15-15-02

HiC-Pro Manual
******************

What is HiC-Pro ?
=================

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to the normalized contact maps. 
The pipeline is flexible, scalable and optimized. It can operate either on a single laptop or on a computational cluster using the PBS-Torque scheduler.

If you use HiC-Pro, please cite :

HiC-Pro: An optimized and flexible pipeline for Hi-C processing. *Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E.*. 2015. submitted

How to install it ?
===================

The HiC-Pro pipeline requires the following dependencies :

* The `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ mapper
* Python (>2.7) with *pysam*, *bx*, *numpy*, and *scipy* libraries
* R with the *RColorBrewer* and *ggplot2* packages
* g++ compiler
* Samtools (>0.1.19)

To install HiC-Pro:

.. code-block:: guess

  tar -zxvf HiC-Pro-master.tar.gz
  cd HiC-Pro-master
  make CONFIG_SYS=config-install.txt install

Note that if some of these dependencies are not installed (i.e. not detected in the $PATH), HiC-Pro will try to install them.
You can also edit the *config-install.txt* file and manually defined the paths to dependencies.

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

In order to process the raw data, HiC-Pro requires three annotation files :

1. A BED file of the restriction fragments after digestion. This file depends both of the restriction enzyme and the reference genome. See the `FAQ <../html/FAQ.html>`_ for details about how to generate this file. A few annotation files are provided with the HiC-Pro sources.

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

2. A table file of chromosomes' size.

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

3. The bowtie2 indexes. See `the bowtie2 manual page <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_ for details about how to create such indexes.

How to use it ?
===============

1. Copy and edit the configuration file *'config-hicpro.txt'* in your local folder. The '[' options are optional.

+---------------+-----------------------------------------+
| SET UP SYSTEM AND PBS/TORQUE MODE                       |
+================+========================================+
| N_CPU          | Number of CPU allows fper job          |
+----------------+----------------------------------------+
| LOGFILE        | Name of the main log file              |
+----------------+----------------------------------------+
| [PBS_SUFFIX]   | Name of PBS/Torque job on the custer   |
+----------------+----------------------------------------+
| [PBS_MEM]      | Memory (RAM) required per job          |
+----------------+----------------------------------------+
| [PBS_WALLTIME] | WallTime allows per job                |
+----------------+----------------------------------------+
| [PBS_MAIL]     | User mail for PBS/Torque report        |
+----------------+----------------------------------------+

------------

+------------------------+---------------------------------------------------------------------------------------------------------------------+
| READS ALIGNMENT OPTIONS                                                                                                                      |
+========================+=====================================================================================================================+
| RAW_DIR                | Link to rawdata folder. The user usually not need to change this option. *Default: rawdata*                         |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| PAIR1_EXT              | Keyword for first mate detection. *Default:_R1*                                                                     |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| PAIR2_EXT              | Keywoard for seconde mate detection. *Default:_R2*                                                                  |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| FORMAT                 | Sequencing qualities encoding. *Default: phred33*                                                                   |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| MIN_MAPQ               | Minimum mapping quality. Reads with lower quality are discarded. *Default: 0*                                       |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| LIGATION SITE          | Ligation site sequence used for reads trimming. Depends on the fill in strategy. *Default: AAGCTAGCTT*              |
+------------------------+---------------------------------------------------------------------------------------------------------------------+ 
| BOWTIE2_IDX_PATH       | Path to bowtie2 indexes                                                                                             |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| BOWTIE2_GLOBAL_OPTIONS | bowtie2 options for mapping step1. *Default: --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder* |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| BOWTIE2_LOCAL_OPTIONS  | bowtie2 options for mapping step2. *Default: --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder* |
+------------------------+---------------------------------------------------------------------------------------------------------------------+

------------

+-----------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| ANNOTATION FILES                                                                                                                                                    |
+=================+===================================================================================================================================================+
| REFERENCE_GENOME| Reference genome prefix used for genome indexes. *Default: hg19*                                                                                  |
+-----------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| GENOME_FRAGMENT | BED file with restriction fragments. Loaded from the ANNOTATION folder in the HiC-Pro installation directory. *Default: HindIII_resfrag_hg19.bed* |
+-----------------+---------------------------------------------------------------------------------------------------------------------------------------------------+
| GENOME_SIZE     | Chromsome size file. Loaded from the ANNOTATION folder in the HiC-Pro installation directory. *Default: chrom_hg19.sizes*                         |
+-----------------+---------------------------------------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| Hi-C PROCESSING                                                                                                                                       |
+=============================+=========================================================================================================================+
| [MIN_INSERT_SIZE]           | Minimum sequenced insert size. Shorter 3C products are discarded                                                        |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| [MAX_INSERT_SIZE]           | Maximum sequenced insert size. Larger 3C products are discarded                                                         |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| GET_ALL_INTERACTION_CLASSES | Create output files with all classes of 3C products. *Default: 1*                                                       |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| GET_PROCESS_BAM             | Create a BAM file with all aligned reads flagged according to their classifaction and mapping category. *Default: 0*    |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_SINGLETON                | Remove singleton reads. *Default: 1*                                                                                    |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_MULTI                    | Remove multi-mapped reads. *Default: 1*                                                                                 |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_DUP                      | Remove duplicated reads' pairs. *Default: 1*                                                                            |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| BIN_SIZE                    | Resolution of contact maps to generate (space separated). *Default: 20000 40000 150000 500000 1000000*                  |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| BIN_STEP                    | Binning step size in ‘n’ coverage _i.e._ window step. *Default: 1*                                                      |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| MATRIX_FORMAT               | Output matrix format. Must be complete, asis, upper or lower. *Default: upper*                                          |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| MAX_ITER                    | Maximum number of iteration for ICE normalization. *Default: 100*                                                       |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| SPARSE_FILTERING            | Define which pourcentage of bins with high sparsity should be force to zero. *Default: 0.02*                            |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| EPS                         | The relative increment in the results before declaring convergence. *Default: 0.1*                                      |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+

------------                                                                                                                                                              

2. Put all fastq files in a rawdata folder. Each fastq file has to be put in a folder per sample.
2. Put all fastq files in a rawdata folder. Each fastq file has to be put in a folder per sample.

3. Run HiC-Pro
   
* Run the complete workflow in a single command line

.. code-block::

   MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE
  
* Run the complete workflow with PBS-Torque

.. code-block:: 

   	MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -p


You will get the following message in the output directory:

.. code-block:: 

  	Please run HiC-Pro in two steps :
  	1- The following command will launch the parallel workflow through 12 torque jobs:
  	qsub HiCPro_step1.sh
  	2- The second command will merge all outputs to generate the contact maps:
  	qsub HiCPro_step2.sh


Execute the displayed command:

.. code-block:: 

  	qsub HiCPro_step1.sh


Then wait for the torque mails... :)
Once executed succesfully (may take several hours), then type:

.. code-block:: 

  	qsub HiCPro_step2.sh


* Run HiC-Pro in sequential mode

HiC-Pro can be run in a step-by-step mode.
Available steps are described in the help command

.. code-block::

  HiC-Pro --help
  usage : HiC-Pro -i INPUT -o OUTPUT -c CONFIG [-s ANALYSIS_STEP] [-p] [-h] [-v]
  Use option -h|--help for more information

  HiC-Pro 2.5.2
  ---------------
  OPTIONS

   -i|--input INPUT : input data folder; Must contains a folder per sample with fastq (or bam) files
   -o|--output OUTPUT : output folder
   -c|--conf CONFIG : configuration file for Hi-C processing
   [-p|--parallel] : if specified run HiC-Pro in PBS/Torque mode
   [-s|--step ANALYSIS_STEP] : run only a subset of the HiC-Pro workflow; if not specified the complete workflow is run
      mapping: perform reads alignment
      proc_hic: perform Hi-C filtering
      quality_checks: run Hi-C quality control plots
      build_contact_maps: build raw inter/intrachromosomal contact maps
      ice_norm : run ICE normalization on contact maps
   [-h|--help]: help
   [-v|--version]: version


As an exemple, if you want to only want to align the sequencing reads, use :

.. code-block::

    	MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -s mapping -s quality_checks

Note that HiC-Pro can be run from aleardy aligned data. In this case, the raw data path (-i) should point into the BAM files.
See te `user"s cases<USER_CASES.rst>`_ for more information


How does HiC-Pro work ?
=======================

.. figure:: images/hicpro_wkflow.png
   :scale: 80%


1. Reads Mapping

Each mate is independantly aligned on the reference genome. The mapping is performed in two steps. First, the reads are aligned using an end-to-end aligner. Second, reads spanning the ligation junction are trimmmed from their 3' end, and aligned on the genome. Aligned reads for both fragment mates are then paired in a single paired-end BAM file. Singletons and multi-hits can be discarded according the confirguration parameters.

2. Fragment assignment and filtering

Each aligned reads can be assigned to one restriction fragment according to the reference genome and the restriction enzyme.
The next step is to separate the invalid ligation products from the valid pairs. Dangling end and self circles pairs are therefore excluded.
Only valid pairs involving two different restriction fragments are used to build the contact maps. Duplicated valid pairs associated to PCR artefacts are discarded.
The fragment assignment can be visualized through a BAM files of aliged pairs where each pair is flagged according to its classification.

3. Quality Controls

HiC-Pro performs a couple of quality controls for most of the analysis steps. The alignment statistics are the first quality controls. Aligned reads in the first (end-to-end) step, and alignment after trimming are reported. Note that in pratice, we ususally observed around 10-20% of trimmed reads. An abnormal level of trimmed reads can reflect a ligation issue.
Once the reads are aligned on the genome, HiC-pro checks the number of singleton, multiple hits or duplicates. The fraction of valid pairs are presented for each type of ligation products. Invalid pairs such as dangling and or self-circle are also represented. A high level of dangling ends, or an imbalance in valid pairs ligation type can be due to a ligation, fill-in or digestion issue.
Finally HiC-Pro also calculated the distribution of fragment size on a subset of valid pairs. Additional statistics will report the fraction of intra/inter-chromosomal contacts, as well as the proportion of short range (<20kb) versus long range (>20kb) contacts.

4. Map builder

Intra et inter-chromosomal contact maps are build for all specified resolutions. The genome is splitted into bins of equal size. Each valid interaction is associated with the genomic bins to generate the raw maps.

5. ICE normalization

Hi-C data can contain several sources of biases which has to be corrected. HiC-Pro proposes a fast implementation of the original ICE normalization algorithm (Imakaev et al. 2012), making the assumption of equal visibility of each fragment. The ICE normalization can be used as a standalone python package and is available `<https://github.com/hiclib/>`_


Data format
===========

A contact map is defined by :

* A list of genomic intervals related to the specified resolution (BED format).
* A matrix, stored as standard triplet sparse format (i.e. list format). Based on the observation that a contact map is symmetric and usually sparse, only non-zero values are stored for half of the matrix. The user can specified if the *'upper'*, *'lower'* or *'complete'* matrix has to be stored. The *'asis'* option allows to store the contacts as they are observed from the valid pairs files.

::

   A   B   10
   A   C   23
   B   C   24
   (...)


This format is memory efficient, and is compatible with other analysis softwares such as the `HiTC Bioconductor package <http://bioconductor.org/packages/release/bioc/html/HiTC.html>`_.







