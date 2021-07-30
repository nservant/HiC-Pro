
# HiC-Pro

### An optimized and flexible pipeline for Hi-C data processing

[![Build Status](https://travis-ci.com/nservant/HiC-Pro.svg?branch=devel_py3)](https://travis-ci.com/nservant/HiC-Pro)

![Conda](https://img.shields.io/badge/Conda-build-brightgreen.svg)
![Singularity](https://img.shields.io/badge/Singularity-build-brightgreen.svg) 
[![Docker](https://img.shields.io/badge/Docker-manual-yellow.svg)](https://hub.docker.com/repository/docker/nservant/hicpro)

![MultiQC](https://img.shields.io/badge/MultiQC-1.8-blue.svg)
[![Forum](https://img.shields.io/badge/Groups-%20join%20chat%20%E2%86%92-4fb99a.svg?style=flat-square)](https://groups.google.com/forum/#!forum/hic-pro)
[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs13059--015--0831--x-lightgrey.svg?style=flat-square)](https://doi.org/10.1186/s13059-015-0831-x)

----

Find documentation and examples at [http://nservant.github.io/HiC-Pro/](http://nservant.github.io/HiC-Pro/)

For any question about HiC-Pro, please contact nicolas.servant@curie.fr or use the [HiC-Pro forum](https://groups.google.com/forum/#!forum/hic-pro)

## What is HiC-Pro ?

HiC-Pro was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to normalized contact maps. It supports the main Hi-C protocols, including digestion protocols as well as protocols that do not require restriction enzymes such as DNase Hi-C. In practice, HiC-Pro was successfully applied to many data-sets including dilution Hi-C, in situ Hi-C, DNase Hi-C, Micro-C, capture-C, capture Hi-C or HiChip data.  
The pipeline is flexible, scalable and optimized. It can operate either on a single laptop or on a computational cluster. HiC-Pro is sequential and each step of the workflow can be run independantly.  
HiC-Pro includes a fast implementatation of the iterative correction method (see the [iced python package](https://github.com/hiclib/iced) for more information).
Finally, HiC-Pro can use phasing data to build [allele-specific contact maps](doc/AS.md).

If you use HiC-Pro, please cite :

*Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E.* HiC-Pro: An optimized and flexible pipeline for Hi-C processing. Genome Biology 2015, 16:259 [doi:10.1186/s13059-015-0831-x](https://doi.org/10.1186/s13059-015-0831-x)

## Containers

### Using HiC-Pro through `conda`

In order to ease the installation of HiC-Pro dependancies, we provide a `yml` file for conda with all required tools.
In order to build your conda environment, first install [miniconda](https://docs.conda.io/en/latest/miniconda.html) and use :

```
conda env create -f MY_INSTALL_PATH/HiC-Pro/environment.yml -p WHERE_TO_INSTALL_MY_ENV
conda activate WHERE_TO_INSTALL_MY_ENV
```

### Using the HiC-Pro `Docker` image

A docker image is automatically build and available on [Docker Hub](https://hub.docker.com/repository/docker/nservant/hicpro)
To pull a Docker image, simply use :

```
docker pull nservant/hicpro:latest
```

Note that the `tag` may depend on the HiC-Pro version.

You can also build your own image from the root folder using

```
docker build -t hicpro:3.1.0 .
```

### Using HiC-Pro through `Singularity`

HiC-Pro provides a Singularity container to ease its installation process.
A ready-to-use container is available [here](https://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_latest_ubuntu.img).

In order to build you own Singularity image;

1- Install singularity

- Linux : http://singularity.lbl.gov/install-linux
- MAC : http://singularity.lbl.gov/install-mac
- Windows : http://singularity.lbl.gov/install-windows

2- Build the singularity HiC-Pro image using the 'Singularity' file available in the HiC-Pro root directory.

```
sudo singularity build hicpro_latest_ubuntu.img MY_INSTALL_PATH/HiC-Pro/envs/Singularity
```

3- Run HiC-pro

You can then either use HiC-Pro using the 'exec' command ;

```
singularity exec hicpro_latest_ubuntu.img HiC-Pro -h
```

Or directly use HiC-Pro within the Singularity shell

```
singularity shell hicpro_latest_ubuntu.img
HiC-Pro -h
```

## How to install it ?

The HiC-Pro pipeline requires the following dependencies :

- The [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) mapper
- Python (>3.7) with *pysam (>=0.15.4)*, *bx-python(>=0.8.8)*, *numpy(>=1.18.1)*, and *scipy(>=1.4.1)* libraries.  
**Note that the current version no longer supports python 2**
- R with the *RColorBrewer* and *ggplot2 (>2.2.1)* packages
- g++ compiler
- samtools (>1.9)
- Unix sort (**which support -V option**) is required ! For Mac OS user, please install the GNU core utilities !

Note that Bowtie >2.2.2 is strongly recommanded for allele specific analysis.  

To install HiC-Pro, be sure to have the appropriate rights and run :

```
tar -zxvf HiC-Pro-master.tar.gz
cd HiC-Pro-master
## Edit config-install.txt file if necessary
make configure
make install
```

Note that if some of these dependencies are not installed (i.e. not detected in the $PATH), HiC-Pro will try to install them.  
You can also edit the *config-install.txt* file and manually defined the paths to dependencies.


|               | SYSTEM CONFIGURATION                                                          |
|---------------|-------------------------------------------------------------------------------|
| PREFIX        | Path to installation folder                                                   |
| BOWTIE2_PATH  | Full path the bowtie2 installation directory                                  |
| SAMTOOLS_PATH | Full path to the samtools installation directory                              |
| R_PATH        | Full path to the R installation directory                                     |
| PYTHON_PATH   | Full path to the python installation directory                                |
| CLUSTER_SYS   | Scheduler to use for cluster submission. Must be TORQUE, SGE, SLURM or LSF    |


## Annotation Files

In order to process the raw data, HiC-Pro requires three annotation files. Note that the pipeline is provided with some Human and Mouse annotation files.  
**Please be sure that the chromosome names are the same than the ones used in your bowtie indexes !**

- **A BED file** of the restriction fragments after digestion. This file depends both of the restriction enzyme and the reference genome. See the [FAQ](doc/FAQ.md) and the [HiC-Pro utilities](doc/UTILS.md) for details about how to generate this file. A few annotation files are provided with the HiC-Pro sources as examples.

```
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
```

- **A table file** of chromosomes' size. This file can be easily find on the UCSC genome browser. Of note, pay attention to the contigs or scaffolds, and be aware that HiC-pro will generate a map per chromosomes pair. For model organisms such as Human or Mouse, which are well annotated, we usually recommand to remove all scaffolds.  

```
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
```

- **The bowtie2 indexes**. See the [bowtie2 manual page](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for details about how to create such indexes.


## How to use it ?

First have a look at the help message !

```
  HiC-Pro --help
  usage : HiC-Pro -i INPUT -o OUTPUT -c CONFIG [-s ANALYSIS_STEP] [-p] [-h] [-v]
  Use option -h|--help for more information

  HiC-Pro 3.1.0
  ---------------
  OPTIONS

   -i|--input INPUT : input data folder; Must contains a folder per sample with input files
   -o|--output OUTPUT : output folder
   -c|--conf CONFIG : configuration file for Hi-C processing
   [-p|--parallel] : if specified run HiC-Pro on a cluster
   [-s|--step ANALYSIS_STEP] : run only a subset of the HiC-Pro workflow; if not specified the complete workflow is run
      mapping: perform reads alignment - require fast files
      proc_hic: perform Hi-C filtering - require BAM files
      quality_checks: run Hi-C quality control plots
      merge_persample: merge multiple inputs and remove duplicates if specified - require .validPairs files
      build_contact_maps: Build raw inter/intrachromosomal contact maps - require .allValidPairs files
      ice_norm : run ICE normalization on contact maps - require .matrix files
   [-h|--help]: help
   [-v|--version]: version
```

- Copy and edit the configuration file *'config-hicpro.txt'* in your local folder. See the [manual](doc/MANUAL.md) for details about the configuration file

- Put all input files in a rawdata folder. The input files have to be organized with **one folder per sample**, such as;

```
   + PATH_TO_MY_DATA
     + sample1
       ++ file1_R1.fastq.gz
       ++ file1_R2.fastq.gz
       ++ ...
     + sample2
       ++ file1_R1.fastq.gz
       ++ file1_R2.fastq.gz
     *...
```

- Run HiC-Pro on your laptop in standalone model

```
    MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_DATA_FOLDER -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE
```

  - Run HiC-Pro on a cluster (TORQUE/SGE/SLURM/LSF)

```
   MY_INSTALL_PATH/bin/HiC-Pro -i FULL_PATH_TO_DATA_FOLDER -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE -p
```

In the latter case, you will have the following message :

```
  Please run HiC-Pro in two steps :
  1- The following command will launch the parallel workflow through 12 torque jobs:
  qsub HiCPro_step1.sh
  2- The second command will merge all outputs to generate the contact maps:
  qsub HiCPro_step2.sh
```

Execute the displayed command from the output folder:

```
  qsub HiCPro_step1.sh
```

Once executed succesfully (may take several hours), run the step using:

```
  qsub HiCPro_step2.sh
```

## Test Dataset

The test dataset and associated results are available [here](https://zerkalo.curie.fr/partage/HiC-Pro/).
Small fastq files (2M reads) extracted from the Dixon et al. 2012 paper are available for test.

```
 ## Get the data. Will download a test_data folder and a configuration file
 wget https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz && tar -zxvf HiCPro_testdata.tar.gz

 ## Edit the configuration file and set the path to Human bowtie2 indexes

 ## Run HiC-Pro
 time HICPRO_INSTALL_DIR/bin/HiC-Pro -c config_test_latest.txt -i test_data -o hicpro_latest_test

Run HiC-Pro 3.1.0
--------------------------------------------
Thu Mar 19, 12:18:10 (UTC+0100)
Bowtie2 alignment step1 ...
Logs: logs/dixon_2M_2/mapping_step1.log
Logs: logs/dixon_2M/mapping_step1.log

--------------------------------------------
Thu Mar 19, 12:18:57 (UTC+0100)
Bowtie2 alignment step2 ...
Logs: logs/dixon_2M_2/mapping_step2.log
Logs: logs/dixon_2M/mapping_step2.log

--------------------------------------------
Thu Mar 19, 12:19:08 (UTC+0100)
Combine R1/R2 alignment files ...
Logs: logs/dixon_2M_2/mapping_combine.log
Logs: logs/dixon_2M/mapping_combine.log

--------------------------------------------
Thu Mar 19, 12:19:13 (UTC+0100)
Mapping statistics for R1 and R2 tags ...
Logs: logs/dixon_2M_2/mapping_stats.log
Logs: logs/dixon_2M/mapping_stats.log

--------------------------------------------
Thu Mar 19, 12:19:15 (UTC+0100)
Pairing of R1 and R2 tags ...
Logs: logs/dixon_2M_2/mergeSAM.log
Logs: logs/dixon_2M/mergeSAM.log

--------------------------------------------
Thu Mar 19, 12:19:25 (UTC+0100)
Assign alignments to restriction fragments ...
Logs: logs/dixon_2M_2/mapped_2hic_fragments.log
Logs: logs/dixon_2M/mapped_2hic_fragments.log

--------------------------------------------
Thu Mar 19, 12:20:10 (UTC+0100)
Merge chunks from the same sample ...
Logs: logs/dixon_2M/merge_valid_interactions.log
Logs: logs/dixon_2M_2/merge_valid_interactions.log

--------------------------------------------
Thu Mar 19, 12:20:11 (UTC+0100)
Merge stat files per sample ...
Logs: logs/dixon_2M/merge_stats.log
Logs: logs/dixon_2M_2/merge_stats.log

--------------------------------------------
Thu Mar 19, 12:20:11 (UTC+0100)
Run quality checks for all samples ...
Logs: logs/dixon_2M/make_Rplots.log
Logs: logs/dixon_2M_2/make_Rplots.log

--------------------------------------------
Thu Mar 19, 12:20:22 (UTC+0100)
Generate binned matrix files ...
Logs: logs/dixon_2M/build_raw_maps.log
Logs: logs/dixon_2M_2/build_raw_maps.log

--------------------------------------------
Thu Mar 19, 12:20:22 (UTC+0100)
Run ICE Normalization ...
Logs: logs/dixon_2M/ice_500000.log
Logs: logs/dixon_2M/ice_1000000.log
Logs: logs/dixon_2M_2/ice_500000.log
Logs: logs/dixon_2M_2/ice_1000000.log

real	2m15,736s
user	4m3,277s
sys	0m24,423s

```
