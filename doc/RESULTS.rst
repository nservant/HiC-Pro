HiC-Pro Results
===============

Mapping results
---------------

The *bowtie_results* folder contains the results of the reads mapping. The results of first mapping step are available in the *bwt2_glob* folder, and the seconnd step in the *bwt2_loc* folder. Final BAM files, reads pairing, and mapping statistics are available on the *bwt2* folder.

List of valid interaction products
----------------------------------

The *hic_results* folder contains all Hi-C processed data, and is further divided in several sub-folders.
The data folder is used to store the valid interaction products (.validPairs), as well as other statisics files.

Intra and inter-chromosomal contact maps
----------------------------------------

The contact maps are then available in the *matrix* folder.

A contact map is defined by :

* A list of genomic intervals related to the specified resolution (BED format).
* A matrix, stored as standard triplet sparse format (i.e. list format). Based on the observation that a contact map is symmetric and usually sparse, only non-zero values are stored for half of the matrix. The user can specified if the *'upper'*, *'lower'* or *'complete'* matrix has to be stored. The *'asis'* option allows to store the contacts as they are observed from the valid pairs files.

::

   A   B   10 
   A   C   23
   B   C   24
   (...)


This format is memory efficient, and is compatible with other analysis softwares such as the `HiTC Bioconductor package <http://bioconductor.org/packages/release/bioc/html/HiTC.html>`_. The *matrix* folder is organized with *raw* and *iced* contact maps for all resolutions.



