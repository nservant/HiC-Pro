## MAKEFILE FOR Hi-C PROCESSING
## Nicolas Servant

## DO NOT EDIT THE REST OF THIS FILE!!

## special characters
comma := ,
space :=
space +=
slash := |
tick := -
undsc := _

all : init mapping proc_hic merge_persample hic_qc build_raw_maps ice_norm

## paralell
all_sub : init mapping proc_hic

init : configure 

mapping: bowtie_global bowtie_local bowtie_combine mapping_stat

proc_hic : bowtie_pairing mapped_2hic_fragments 

## per sample - no parallelisation
all_persample : merge_persample build_raw_maps hic_qc ice_norm

merge_persample: merge_valid_interactions merge_allstat

quality_checks : hic_qc

build_contact_maps: build_raw_maps


######################################
## System
##
######################################
config_check:
ifndef CONFIG_FILE
	$(error CONFIG_FILE is not defined)
else		
include $(CONFIG_FILE)
endif
ifndef CONFIG_SYS
	$(error CONFIG_SYS is not defined)
else		
include $(CONFIG_SYS)
endif

make_cluster_script: config_check init
ifdef MAKE_OPTS
	@$(CLUSTER_SCRIPT) -c $(CONFIG_FILE) -s $(MAKE_OPTS)
else
	@$(CLUSTER_SCRIPT) -c $(CONFIG_FILE)
endif

clean:
	@/bin/rm -rf $(LOGS_DIR)/*
	@/bin/rm -rf $(BOWTIE2_OUTPUT_DIR)/bwt2_glob  
	@/bin/rm -rf $(BOWTIE2_OUTPUT_DIR)/bwt2_loc  


debug: config_check
	@echo "RAW_DIR="$(RAW_DIR)
	@echo "FASTQ_FILE="$(READSFILE_FQ)

######################################
## Configure outputs
##
######################################

## Create output folders
configure:  config_check
	@mkdir -p $(TMP_DIR)
	@mkdir -p $(LOGS_DIR)

######################################
## Bowtie2 Global Alignment
##
######################################

## End-to-end Alignement
bowtie_global:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Bowtie2 alignment step1 ..." 
	@$(SCRIPTS)/bowtie_wrap.sh -c $(CONFIG_FILE) -m step1 -u
	@echo

######################################
##  Bowtie2 Local Alignment
##
######################################

## Cut-and-map Alignement
bowtie_local:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Bowtie2 alignment step2 ..."
	@$(SCRIPTS)/bowtie_wrap.sh -c $(CONFIG_FILE) -m step2
	@echo


######################################
## Merge Bowtie2 local and global mapping
## 
######################################

## Merge global and local alignment in a single file
bowtie_combine:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Combine R1/R2 alignment files ..."
	@$(SCRIPTS)/bowtie_combine.sh -c $(CONFIG_FILE) 
	@echo

## Compute mapping statistics
mapping_stat:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Mapping statistics for R1 and R2 tags ..." 
	@$(SCRIPTS)/mapping_stat.sh -c $(CONFIG_FILE)
	@echo

mapping_clean: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Cleaning mapping files ..."
	@rm -rf $(BOWTIE2_OUTPUT_DIR)/bwt2_global
	@rm -rf $(BOWTIE2_OUTPUT_DIR)/bwt2_local
	@echo

######################################
## Hi-C processing 
##
######################################

## Pairing of R1 and R2 mates and reads filtering
bowtie_pairing:  config_check
	@echo "--------------------------------------------"
	@date 
	@echo "Pairing of R1 and R2 tags ..." 
	@$(SCRIPTS)/bowtie_pairing.sh -c $(CONFIG_FILE)
	@echo

## Assign alignments to restriction fragments
mapped_2hic_fragments:  config_check
	@echo "--------------------------------------------"
	@date 
	@echo "Assign alignments to restriction fragments ..." 
	@$(SCRIPTS)/mapped_2hic_fragments.sh -c $(CONFIG_FILE)
	@echo

## Combine multiple BAM files from the same sample, and remove duplicates
merge_valid_interactions:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Merge chunks from the same sample ..." 
	@$(SCRIPTS)/merge_valid_interactions.sh -c $(CONFIG_FILE)
	@echo

merge_allstat: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Merge stat files per sample ..." 
	@$(SCRIPTS)/merge_stats.sh -c $(CONFIG_FILE)
	@echo

build_raw_maps:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Generate binned matrix files ..." 
	@$(SCRIPTS)/build_raw_maps.sh -c $(CONFIG_FILE)
	@echo

hic_qc:
	@echo "--------------------------------------------" 
	@date 
	@echo "Run quality checks for all samples ..."
	@$(SCRIPTS)/make_plots.sh -c $(CONFIG_FILE) -p "all"
	@echo

######################################
## Normalization
##
######################################

## Apply ICE normalization
ice_norm:
	@echo "--------------------------------------------" 
	@date 
	@echo "Run ICE Normalization ..." 
	@$(SCRIPTS)/ice_norm.sh -c $(CONFIG_FILE)

