.. _MANUAL:

.. Nicolas Servant
.. HiC-Pro
.. v2.7.0
.. 15-09-02

HiC-Pro Manual
******************
Modified - 10th September 2017
Reference version - HiC-Pro 2.8.0

Setting the configuration file
==============================

1. Copy and edit the configuration file *'config-hicpro.txt'* in your local folder. The '[]' options are optional and can be undefined.

+----------------+----------------------------------------+
| SET UP SYSTEM AND CLUSTER MODE                          |
+================+========================================+
| N_CPU          | Number of CPU allows per job           |
+----------------+----------------------------------------+
| LOGFILE        | Name of the main log file              |
+----------------+----------------------------------------+
| [JOB_NAME  ]   | Name of the job on the custer          |
+----------------+----------------------------------------+
| [JOB_MEM]      | Memory (RAM) required per job          |
+----------------+----------------------------------------+
| [JOB_WALLTIME] | WallTime allows per job                |
+----------------+----------------------------------------+
| [JOB_MAIL]     | User mail for PBS/Torque report        |
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
| BOWTIE2_IDX_PATH       | Path to bowtie2 indexes                                                                                             |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| BOWTIE2_GLOBAL_OPTIONS | bowtie2 options for mapping step1. *Default: --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder* |
+------------------------+---------------------------------------------------------------------------------------------------------------------+
| BOWTIE2_LOCAL_OPTIONS  | bowtie2 options for mapping step2. *Default: --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder* |
+------------------------+---------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| ANNOTATION FILES      |                                                                                                                                   |
+=======================+===================================================================================================================================+
| REFERENCE_GENOME      | Reference genome prefix used for genome indexes. *Default: hg19*                                                                  |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| GENOME_SIZE           | Chromsome size file. Loaded from the ANNOTATION folder in the HiC-Pro installation directory. *Default: chrom_hg19.sizes*         |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| [CAPTURE_BED]         | BED file of target regions to focus on (mainly used for capture Hi-C data                                                         |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------+
| [ALLELE_SPECIFIC_SNP] | VCF file to SNPs which can be used to distinguish parental origin. See the :ref:`allele specific section <AS>` for more details   |
+-----------------------+-----------------------------------------------------------------------------------------------------------------------------------+

------------

+---------------------------+----------------------------------------------------------------------------------------------------------------------------------+
| ALLLELE SPECIFIC ANALYSIS |                                                                                                                                  |
+=======================+======================================================================================================================================+
| [ALLELE_SPECIFIC_SNP]     | VCF file to SNPs which can be used to distinguish parental origin. See the :ref:`allele specific section <AS>` for more details  |
+---------------------------+----------------------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------+
| DIGESTION Hi-C        |                                                                                                                                          |
+=======================+==========================================================================================================================================+
| [GENOME_FRAGMENT]     | BED file with restriction fragments. Full path or name of file available in the ANNOTATION folder. *Default: HindIII_resfrag_hg19.bed*   |
+-----------------------+------------------------------------------------------------------------------------------------------------------------------------------+
| [LIGATION SITE]       | Ligation site sequence used for reads trimming. Depends on the fill in strategy. *Example: AAGCTAGCTT*                                   |
+------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+ 
| [MIN_FRAG_SIZE]       | Maximum size of restriction fragments to consider for the Hi-C processing. *Example: 100*                                                | 
+------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+ 
| [MAX_FRAG_SIZE]       | Maximum size of restriction fragments to consider for the Hi-C processing. *Example: 100000*                                             |
+------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+ 
| [MIN_INSERT_SIZE]     | Minimum sequenced insert size. Shorter 3C products are discarded. *Example: 100*                                                         |
+-----------------------------+------------------------------------------------------------------------------------------------------------------------------------+
| [MAX_INSERT_SIZE]     | Maximum sequenced insert size. Larger 3C products are discarded. *Example: 600*                                                          |
+-----------------------------+------------------------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| Hi-C PROCESSING                                                                                                                                       |
+=============================+=========================================================================================================================+
| [MIN_CIS_DIST]              | Filter short range contact below the specified distance. Mainly useful for DNase Hi-C. *Example: 1000*                  |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| GET_ALL_INTERACTION_CLASSES | Create output files with all classes of 3C products. *Default: 0*                                                       |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| GET_PROCESS_BAM             | Create a BAM file with all aligned reads flagged according to their classifaction and mapping category. *Default: 0*    |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_SINGLETON                | Remove singleton reads. *Default: 1*                                                                                    |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_MULTI                    | Remove multi-mapped reads. *Default: 1*                                                                                 |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| RM_DUP                      | Remove duplicated reads' pairs. *Default: 1*                                                                            |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| GENOME-WIDE CONTACT MAPS                                                                                                                              |
+=============================+=========================================================================================================================+
| BIN_SIZE                    | Resolution of contact maps to generate (space separated). *Default: 20000 40000 150000 500000 1000000*                  |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| BIN_STEP                    | Binning step size in ‘n’ coverage _i.e._ window step. *Default: 1*                                                      |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+
| MATRIX_FORMAT               | Output matrix format. Must be complete, upper. *Default: upper*. *Deprecated: asis, lower*                              |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------+

------------

+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+
| NORMALIZATION                                                                                                                                               |
+===================================+=========================================================================================================================+
| MAX_ITER                          | Maximum number of iteration for ICE normalization. *Default: 100*                                                       |
+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+
| SPARSE_FILTERING - **deprecated** | Define which pourcentage of bins with high sparsity should be force to zero. *Default: 0.02*                            |
+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+
| FILTER_LOW_COUNT_PERC             | Define which pourcentage of bins with low counts should be force to zero. *Default: 0.02*. Replace SPARSE_FILTERING     |
+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+
| FILTER_HIGH_COUNT_PERC            | Define which pourcentage of bins with low counts should be discarded before normalization. *Default: 0*                 |
+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+
| EPS                               | The relative increment in the results before declaring convergence. *Default: 0.1*                                      |
+-----------------------------------+-------------------------------------------------------------------------------------------------------------------------+

------------           

                                                                                                                                                   

Run HiC-Pro in sequential mode
==============================

HiC-Pro can be run in a step-by-step mode.
Available steps are described in the help command.

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
   [-p|--parallel] : if specified run HiC-Pro on a cluster
   [-s|--step ANALYSIS_STEP] : run only a subset of the HiC-Pro workflow; if not specified the complete workflow is run
      mapping: perform reads alignment
      proc_hic: perform Hi-C filtering
      quality_checks: run Hi-C quality control plots
      build_contact_maps: build raw inter/intrachromosomal contact maps
      ice_norm: run ICE normalization on contact maps
   [-h|--help]: help
   [-v|--version]: version


As an exemple, if you want to only want to only align the sequencing reads and run a quality control, use :

.. code-block:: guess

    	MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_RAW_DATA -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -s mapping -s quality_checks

Note that in sequential mode, the INPUT argument depends on the analysis step. See te :ref:`user's cases <USERCASES>` for more examples.

+-----------------------+--------------------+
| INPUT DATA TYPE IN STEPWISE MODE           |
+=======================+====================+
|  -s mapping           | .fastq(.gz) files  |
+-----------------------+--------------------+
| -s proc_hic           | .bam files         |
+-----------------------+--------------------+
| -s quality_checks     | .bam files         |
+-----------------------+--------------------+
| -s merge_persample    | .validPairs files  |
+-----------------------+--------------------+
| -s build_contact_maps | .validPairs files  |
+-----------------------+--------------------+
| -s ice_norm           | .matrix files      |
+-----------------------+--------------------+


How does HiC-Pro work ?
=======================

The HiC-Pro workflow can be divided in five main steps presented below.

.. figure:: images/hicpro_wkflow.png
   :scale: 80%


1. **Reads Mapping**

| Each mate is independantly aligned on the reference genome. The mapping is performed in two steps. First, the reads are aligned using an end-to-end aligner. Second, reads spanning the ligation junction are trimmmed from their 3' end, and aligned back on the genome. Aligned reads for both fragment mates are then paired in a single paired-end BAM file. Singletons and multi-hits can be discarded according the confirguration parameters. Note that if if the *LIGATION_SITE* parameter in the not defined, HiC-Pro will skip the second step of mapping.

2. **Fragment assignment and filtering**

| Each aligned reads can be assigned to one restriction fragment according to the reference genome and the restriction enzyme.
| The next step is to separate the invalid ligation products from the valid pairs.
| Here is the list of pairs classified as invalid by HiC-Pro :

* Dangling end, i.e. unligated fragments (both reads mapped on the same restriction fragment)
* Self circles, i.e. fragments ligated on themselves (both reads mapped on the same restriction fragment in inverted orientation
* Religation, i.e. ligation of juxtaposed fragments
* Dumped pairs, i.e. any pairs that do not match the filtering criteria on inserts size, restriction fragments size or for which we were not able to reconstruct the ligation product.

| Only valid pairs involving two different restriction fragments are used to build the contact maps. Duplicated valid pairs associated to PCR artefacts are discarded.
| The fragment assignment can be visualized through a BAM files of aliged pairs where each pair is flagged according to its classification.
| In case of Hi-C protocols that do not require a restriction enzyme such as DNase Hi-C or micro Hi-C, the assignment to a restriction is not possible. If no *GENOME_FRAGMENT* file are specified, this step is ignored. Short range interactions can however still be discarded using the *MIN_CIS_DIST* parameter.

3. **Quality Controls**

| HiC-Pro performs a couple of quality controls for most of the analysis steps. The alignment statistics are the first quality controls. Aligned reads in the first (end-to-end) step, and alignment after trimming are reported. Note that in pratice, we ususally observed around 10-20% of trimmed reads. An abnormal level of trimmed reads can reflect a ligation issue.
| Once the reads are aligned on the genome, HiC-pro checks the number of singleton, multiple hits or duplicates. The fraction of valid pairs are presented for each type of ligation products. Invalid pairs such as dangling and or self-circle are also represented. A high level of dangling ends, or an imbalance in valid pairs ligation type can be due to a ligation, fill-in or digestion issue.
| Finally HiC-Pro also calculated the distribution of fragment size on a subset of valid pairs. Additional statistics will report the fraction of intra/inter-chromosomal contacts, as well as the proportion of short range (<20kb) versus long range (>20kb) contacts.

4. **Map builder**

| Intra et inter-chromosomal contact maps are build for all specified resolutions. The genome is splitted into bins of equal size. Each valid interaction is associated with the genomic bins to generate the raw maps.

5. **ICE normalization**

| Hi-C data can contain several sources of biases which has to be corrected. HiC-Pro proposes a fast implementation of the original ICE normalization algorithm (Imakaev et al. 2012), making the assumption of equal visibility of each fragment. The ICE normalization can be used as a standalone python package through the `iced python package <https://github.com/hiclib/>`_.


Browsing the results
====================

All outputs follow the input organization, with one folder per sample.
See the :ref:`results <RES>` section for more information.

* *bowtie_results*

The *bowtie_results* folder contains the results of the reads mapping. The results of first mapping step are available in the *bwt2_glob* folder, and the seconnd step in the *bwt2_loc* folder. Final BAM files, reads pairing, and mapping statistics are available on the *bwt2* folder. Note that once HiC-Pro has been run, all files in *bwt2_glob* or *bwt2_loc* folders can be removed. These files take a significant amount of disk space and are not useful anymore.

* *hic_results*

| This folder contains all Hi-C processed data, and is further divided in several sub-folders.
| The *data* folder is used to store the valid interaction products (*.validPairs*), as well as other statisics files.

| The *validPairs* are stored using a simple tab-delimited text format ;
| read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 / strand_reads2 / fragment_size [/ allele_specific_tag]
| One *validPairs* file is generated per reads chunck. These files are then merged in the *allValidPairs*, and duplicates are removed if specified in the configuration file.

| The contact maps are then available in the *matrix* folder. The *matrix* folder is organized with *raw* and *iced* contact maps for all resolutions.
| Contact maps are stored as a triplet sparse format ;
| bin_i / bin_j / counts_ij
| Only no zero values are stored. BED file described the genomic bins are also generated. Note that *abs* and *ord* files are identical in the context of Hi-C data as the contact maps are symmetric.

| Finally, the *pic* folder contains graphical outputs of the quality control checks.






