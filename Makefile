## MAKEFILE FOR Hi-C PROCESSING
## Nicolas Servant

## DO NOT EDIT THE REST OF THIS FILE!!
include $(CONFIG_FILE)

## special characters
comma := ,
space :=
space +=
slash := |
tick := -
undsc := _

#####################################
## Define Variables
##
#####################################

ALL_CHRS := $(subst $(space),$(comma),$(patsubst %,chr%,$(MOUSE_CHRMS)))

READSFILE_FQ := $(wildcard $(RAW_DIR)/*.f*q)
READSFILE_FQ += $(wildcard $(RAW_DIR)/*.f*q.gz)
READSFILE_FQ_R1 := $(wildcard $(RAW_DIR)/*_R1_*.f*q)
READSFILE_FQ_R1 += $(wildcard $(RAW_DIR)/*_R1_*.f*q.gz)
READSFILE_FQ_R2 := $(wildcard $(RAW_DIR)/*_R2_*.f*q)
READSFILE_FQ_R2 += $(wildcard $(RAW_DIR)/*_R2_*.f*q.gz)

RES_FILE_NAME := $(notdir $(RAW_DIR))
RES_FILE_NAME_OBJ := $(subst $(tick),$(undsc),$(RES_FILE_NAME))

BOWTIE_GENOME_LOGFILE := $(LOGS_DIR)/bowtie_mapping

BOWTIE2_IDX = $(BOWTIE2_IDX_PATH)/$(ORGANISM)
BOWTIE2_LOCAL_OUTPUT_DIR = $(BOWTIE2_OUTPUT_DIR)/bwt2_local/$(RES_FILE_NAME)
BOWTIE2_GLOBAL_OUTPUT_DIR = $(BOWTIE2_OUTPUT_DIR)/bwt2_global/$(RES_FILE_NAME)
BOWTIE2_FINAL_OUTPUT_DIR = $(BOWTIE2_OUTPUT_DIR)/bwt2/$(RES_FILE_NAME)
UNMAP_READ_DIR = $(BOWTIE2_OUTPUT_DIR)/unmap/$(RES_FILE_NAME)
DATA_DIR = $(MAPC_OUTPUT)/data/$(RES_FILE_NAME)
##DOC_DIR = $(MAPC_OUTPUT)/doc/$(RES_FILE_NAME)
PIC_DIR = $(MAPC_OUTPUT)/pic/$(RES_FILE_NAME)
##PROB_DIR = $(MAPC_OUTPUT)/prob_matrix/$(RES_FILE_NAME)
MATRIX_DIR = $(MAPC_OUTPUT)/matrix/$(RES_FILE_NAME)/raw
ICED_MATRIX_DIR = $(MAPC_OUTPUT)/matrix/$(RES_FILE_NAME)/iced
LGF_MATRIX_DIR = $(MAPC_OUTPUT)/matrix/$(RES_FILE_NAME)/lgf
RDATA_DIR = $(MAPC_OUTPUT)/rdata/$(RES_FILE_NAME)


GLOBAL_SAM=$(patsubst %, $(BOWTIE2_GLOBAL_OUTPUT_DIR)/%_$(ORGANISM).bwt2glob.sam, $(basename $(notdir $(READSFILE_FQ))))
GLOBAL_BAM=$(GLOBAL_SAM:.sam=.bam)
GLOBAL_PROCESS=$(GLOBAL_SAM:.sam=.aln)

PROCESS_PREFIX=$(patsubst %, $(BOWTIE2_GLOBAL_OUTPUT_DIR)/%_$(ORGANISM).bwt2glob.sam, $(basename $(notdir $(READSFILE_FQ))))


all : init mapping proc_hic build_contact_maps clean

init : configure src_compile

mapping: bowtie_global bowtie_local merge_global_local mapping_stat ##plot_MappingProportion 

proc_hic : bowtie_pairing mapped2HiCFragments merge_rmdup

build_contact_maps: build_raw_maps ##matrix2RData


##mapping_prep_hic: mapping_noplot prep_hic

##mapping_noplot: mapping_glob mapping_loc merge_global_local

##prep_hic: separateAlignment_all Mapped2HiCFragments_all

##mapping_glob: bowtie_global bowtie_global_stat bowtie_global_pp

##mapping_loc:  bowtie_local bowtie_local_stat bowtie_local_pp 

##norm: ICEnorm LGFnorm

test_config:
	$(SCRIPTS)/test_config.sh $(CONFIG_FILE)

debug:
	@echo $(RAW_DIR)
	@echo $(READSFILE_FQ)


######################################
## System
##
######################################

make_torque_script:
	@$(SCRIPTS)/make_torque_script.sh -c $(CONFIG_FILE) $(TORQUE_SUFFIX)

clean:
	/bin/rm -f $(BOWTIE2_GLOBAL_OUTPUT_DIR)/*.sam
	/bin/rm -f $(BOWTIE2_LOCAL_OUTPUT_DIR)/*.sam

reset: 
ifdef LOGS_DIR
	/bin/rm -rf bowtie_results hic_results $(LOGS_DIR)/*
endif

######################################
## Configure outputs
##
######################################

configure:
ifndef CONFIG_FILE
	$(error CONFIG_FILE is not defined)
endif
	mkdir -p $(BOWTIE2_OUTPUT_DIR)
	mkdir -p $(MAPC_OUTPUT)
	mkdir -p $(TMP_DIR)
	mkdir -p $(DATA_DIR)
	mkdir -p $(PIC_DIR)
	mkdir -p $(MATRIX_DIR)
	mkdir -p $(ICED_MATRIX_DIR)
	mkdir -p $(LOGS_DIR)
	@echo "## Hi-C Mapping $(VERSION)" > $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "\n## INPUT FILES :" >> $(LOGFILE)
	@echo $(READSFILE_FQ)  >> $(LOGFILE)

######################################
## Compile
##
######################################

src_compile: $(SOURCES)/build_matrix.cpp
	(cd $(SOURCES); g++ -Wall -O2 -std=c++0x -o build_matrix build_matrix.cpp)


######################################
## Bowtie2 Global Alignment
##
######################################

## Global Alignement
bowtie_global:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 global alignment ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_wrap.sh -c $(CONFIG_FILE) -u

## Bowtie post-processing
# bowtie_global_pp:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)	
# 	@echo "Bowtie2 global alignment post-process ..." >> $(LOGFILE)
# 	$(SCRIPTS)/mapping_pp.sh -c $(CONFIG_FILE) -u

## bowtie2 global mapping stats
# bowtie_global_stat:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Bowtie2 global mapping statistics for R1 and R2 tags ..." >> $(LOGFILE)
# 	$(SCRIPTS)/mapping_stat.sh -c $(CONFIG_FILE)

######################################
##  Bowtie2 Local Alignment
##
######################################

## Local Alignement
bowtie_local:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 local alignment ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_wrap.sh -c $(CONFIG_FILE) -l

## Bowtie post-processing
# bowtie_local_pp:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)	
# 	@echo "Bowtie2 local alignment post-process ..." >> $(LOGFILE)
# 	$(SCRIPTS)/mapping_pp.sh -c $(CONFIG_FILE) -l

## bowtie2 global mapping stats
# bowtie_local_stat:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Bowtie2 local mapping statistics for R1 and R2 tags ..." >> $(LOGFILE)
# 	$(SCRIPTS)/mapping_stat.sh -c $(CONFIG_FILE) -l

######################################
## Merge Bowtie2 local and global mapping
## 
######################################

merge_global_local:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Merge both alignment ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_combine.sh -c $(CONFIG_FILE)

mapping_stat:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 mapping statistics for R1 and R2 tags ..." >> $(LOGFILE)
	$(SCRIPTS)/mapping_stat.sh -c $(CONFIG_FILE)

# plot_MappingProportion:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Plot Mapping Proportion ..." >> $(LOGFILE)
# 	$(SCRIPTS)/plotMappingPortion.sh -c $(CURDIR)/$(CONFIG_FILE)


######################################
## Hi-C processing 
##
######################################

## Pairing of R1 and R2 mates and reads filtering
bowtie_pairing:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Pairing of R1 and R2 tags ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_pairing.sh -c $(CONFIG_FILE)

## Assign alignments to regions segmented by HindIII sites
mapped2HiCFragments:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Assign alignments to HindIII sites ..." >> $(LOGFILE)
	$(SCRIPTS)/overlapMapped2HiCFragments.sh -c $(CONFIG_FILE) > $(LOGS_DIR)/overlapRS.log

merge_rmdup:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Merge multiple files from the same sample ..." >> $(LOGFILE)
	$(SCRIPTS)/mergeValidInteractions.sh -c $(CONFIG_FILE) > $(LOGS_DIR)/mergeMulti.log

merge_rmdup2:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Merge multiple files from the same sample ..." >> $(LOGFILE)
	$(SCRIPTS)/mergeValidInteractions_v2.sh -c $(CONFIG_FILE) > $(LOGS_DIR)/mergeMulti_v2.log

## Only keep alignments which have both ends uniquely mapped to the genome
## Inputs : final.aln files (R1 + R2)
## Outputs : out files
#separateAlignment_all:
#	@echo "--------------------------------------------" >> $(LOGFILE)
#	@date >> $(LOGFILE)
#	@echo "Filter out reads alignement ..." >> $(LOGFILE)
#	 $(SCRIPTS)/separatePEalignment.sh -c $(CONFIG_FILE)

## Assign alignments to regions segmented by HindIII sites
# Mapped2HiCFragments_all:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Assign alignments to HindIII sites ..." >> $(LOGFILE)
# 	$(SCRIPTS)/overlapMapped2HiCFragments.sh -c $(CONFIG_FILE) # > $(LOGS_DIR)/overlapRS.log

## assign segmented regions to defined bins
build_raw_maps:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Generate binned matrix files ..." >> $(LOGFILE)
	$(SCRIPTS)/assignRead2bins.sh -c $(CONFIG_FILE) > $(LOGS_DIR)/build_raw_maps.log ##-i $(DATA_DIR)/$(RES_FILE_NAME).$(ORGANISM).interaction -g $(ANNOT_DIR)/chrom.sizes

# plot_FragmentInfo:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Plot Fragment Informations ..." >> $(LOGFILE)
# 	$(SCRIPTS)/plotHiCFragment.sh -c $(CURDIR)/$(CONFIG_FILE)

# matrix2RData:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Generate RData files ..." >> $(LOGFILE)
# 	$(SCRIPTS)/matrix2RData.sh -c $(CURDIR)/$(CONFIG_FILE)

######################################
## Normalization
##
######################################

# ICEnorm:
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Run ICE Normalization ..." >> $(LOGFILE)
# 	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),$(R_PATH)/R --no-save CMD BATCH "--args rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE).RData' norm='ICE' cpu='$(N_CPU)' outDir='$(RDATA_DIR)'" $(SCRIPTS)/hicNorm.R $(LOGS_DIR)/hicNorm.Rout;)
# 	@echo "Generate Matrix files from RData ..." >> $(LOGFILE)
# 	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),mkdir -p $(ICED_MATRIX_DIR)/$(BSIZE); $(R_PATH)/R --no-save CMD BATCH "--args matDir='$(ICED_MATRIX_DIR)/$(BSIZE)' cpu='$(N_CPU)' rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_iced.RData' org='$(ORGANISM)'" $(SCRIPTS)/rData2matrix.R $(LOGS_DIR)/rData2matrix.Rout;)


# LGFnorm: 
# 	@echo "--------------------------------------------" >> $(LOGFILE)
# 	@date >> $(LOGFILE)
# 	@echo "Run LGF Normalization ..." >> $(LOGFILE)
# 	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),$(R_PATH)/R --no-save CMD BATCH "--args rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_annot.RData' norm='LGF' cpu='$(N_CPU)' outDir='$(RDATA_DIR)'" $(SCRIPTS)/hicNorm.R $(LOGS_DIR)/hicNorm.Rout;)
# 	@echo "Generate Matrix files from RData ..." >> $(LOGFILE)
# 	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),mkdir -p $(LGF_MATRIX_DIR)/$(BSIZE); $(R_PATH)/R --no-save CMD BATCH "--args matDir='$(LGF_MATRIX_DIR)/$(BSIZE)' cpu='$(N_CPU)' rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_annot_lgf.RData' org='$(ORGANISM)'" $(SCRIPTS)/rData2matrix.R $(LOGS_DIR)/rData2matrix.Rout;)

