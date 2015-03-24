HiC-Pro Utilities
=================

..  _UTILS:


HiC-Pro provides a couple of utilities which are not part of the main pipeline, but was designed to help the user to perform a couple of tasks.
All utilities are installed in /bin/utils/.

**1- split_reads.py - How can I split my reads in chunck ?**

Split reads is highly recommanded in HiC-pro parallel mode.
To split the reads in chunck, you can simply used the split_reads.py utility.

.. code-block:: bash

   ## split the fastq file in blocks of 10M reads
   HICPRO_PATH/bin/utils/split_reads.py --results_folder OUTPUT --nreads READS_NB INPUT_FASTQ




