User cases
==========

Let's defined a few variables ;
* TEST_DATA is my input folder with one folder per sample and all my fastq files
* RES_PREFIX is the prefix of my output directory
* config_test.txt is my configuration file

1. Running HiC-Pro's complete pipeline on my laptop

.. code-block:: guess

   HiC-Pro -i ${TEST_DATA} -o ${RES_PREFIX}_1 -c config_test.txt


2. Running HiC-Pro on a cluster with PBS

.. code-block:: guess

   HiC-Pro -i ${TEST_DATA} -o ${RES_PREFIX}_2 -c config_test.txt -p
   cd ${TEST_DATA}
   qsub HiCPro_step1_test.sh
   ## wait ...
   qsub HiCPro_step2_test.sh


3. Running HiC-Pro using stepwise function

I first want to align to raw reads and make the quality controls on the mapping results using the parallel mode

.. code-block:: guess

   HiC-Pro -i ${TEST_DATA} -o ${RES_PREFIX}_3 -c config_test.txt -s mapping -s quality_checks -p
   cd ${RES_PREFIX}_3
   qsub HiCPro_step1_mappingtest.sh


Once the reads are mapped, let's run the Hi-C processing step from the aligned data on my own laptop

.. code-block:: guess

   HiC-Pro -i ${RES_PREFIX}_3/bowtie_results/bwt2 -o ${RES_PREFIX}_3.1 -c config_test.txt -s proc_hic -s quality_checks


Everything looks good, I can then create the contact maps and normalize them

.. code-block:: guess

   HiC-Pro -i ${RES_PREFIX}_3/bowtie_results/bwt2 -o ${RES_PREFIX}_3.2 -c config_test.txt -s proc_hic -s build_contact_maps -s ice_norm


4. To run an allele specific analysis, first set the ALLELE_SPECIFIC_SNP variable in the configuration file and check whether the N-masked genome is available. Then simply run HiC-Pro as usually.

.. code-block:: guess

   HiC-Pro -i ${TEST_DATA} -o ${RES_PREFIX}_5 -c config_test_as.txt


