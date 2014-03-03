## MAKEFILE FOR Hi-C PROCESSING
## Chong-jian Chen
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
DOC_DIR = $(MAPC_OUTPUT)/doc/$(RES_FILE_NAME)
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


all : configure mapping hic norm

mapping: mapping_glob mapping_loc merge_global_local plot_MappingProportion mapping_clean

hic: separateAlignment_all Mapped2HiCFragments_all plot_FragmentInfo AssignRead2bin_all  matrix2RData

mapping_glob: bowtie_global bowtie_global_pp bowtie_global_stat

mapping_loc:  bowtie_local bowtie_local_pp bowtie_local_stat

norm: ICEnorm LGFnorm

debug:
	@echo $(RAW_DIR)
	@echo $(READSFILE_FQ)

######################################
## Configure outputs
##
######################################

configure:
	mkdir -p $(BOWTIE2_OUTPUT_DIR)
	mkdir -p $(BOWTIE2_GLOBAL_OUTPUT_DIR)
	mkdir -p $(BOWTIE2_LOCAL_OUTPUT_DIR)
	mkdir -p $(BOWTIE2_FINAL_OUTPUT_DIR)
	mkdir -p $(MAPC_OUTPUT)
	mkdir -p $(DATA_DIR)
	mkdir -p $(DOC_DIR)
	mkdir -p $(PIC_DIR)
	mkdir -p $(MATRIX_DIR)
	mkdir -p $(ICED_MATRIX_DIR)
	mkdir -p $(RDATA_DIR)
	mkdir -p $(LOGS_DIR)
	@echo "## Hi-C Mapping v1.0" > $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "\n## INPUT FILES :" >> $(LOGFILE)
	@echo $(READSFILE_FQ)  >> $(LOGFILE)

######################################
## Bowtie2 Global Alignment
##
######################################

## Global Alignement
bowtie_global:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 global alignment ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_wrap.sh -c $(CONFIG_FILE)

## Bowtie post-processing
bowtie_global_pp:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)	
	@echo "Bowtie2 global alignment post-process ..." >> $(LOGFILE)
	$(SCRIPTS)/mapping_pp.sh -i $(BOWTIE2_GLOBAL_OUTPUT_DIR) -u


## bowtie2 global mapping stats
bowtie_global_stat:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 global mapping statistics for R1 and R2 tags ..." >> $(LOGFILE)
	$(SCRIPTS)/mappingstat.sh -i $(BOWTIE2_GLOBAL_OUTPUT_DIR) -o ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${RES_FILE_NAME}

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
bowtie_local_pp:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)	
	@echo "Bowtie2 local alignment post-process ..." >> $(LOGFILE)
	$(SCRIPTS)/mapping_pp.sh -i $(BOWTIE2_LOCAL_OUTPUT_DIR)

## bowtie2 global mapping stats
bowtie_local_stat:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Bowtie2 local mapping statistics for R1 and R2 tags ..." >> $(LOGFILE)
	$(SCRIPTS)/mappingstat.sh -i $(BOWTIE2_LOCAL_OUTPUT_DIR) -o ${BOWTIE2_LOCAL_OUTPUT_DIR}/${RES_FILE_NAME}

######################################
## Merge Bowtie2 local and global mapping
## 
######################################

merge_global_local:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Merge both alignment ..." >> $(LOGFILE)
	$(SCRIPTS)/bowtie_combine.sh -c $(CONFIG_FILE) -g $(BOWTIE2_GLOBAL_OUTPUT_DIR) -l $(BOWTIE2_LOCAL_OUTPUT_DIR) -o $(BOWTIE2_FINAL_OUTPUT_DIR)

plot_MappingProportion:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Plot Mapping Proportion ..." >> $(LOGFILE)
	$(R_PATH)/R --no-save CMD BATCH "--args picDir='$(PIC_DIR)' docDir='$(DOC_DIR)' bwtDir='$(BOWTIE2_OUTPUT_DIR)' sampleName='$(RES_FILE_NAME)'" $(SCRIPTS)/plotMappingPortion.R $(LOGS_DIR)/plotMappingPortion.Rout

mapping_clean:
	/bin/rm -f $(BOWTIE2_GLOBAL_OUTPUT_DIR)/*.sam
	/bin/rm -f $(BOWTIE2_LOCAL_OUTPUT_DIR)/*.sam
	/bin/rm -f $(BOWTIE2_FINAL_OUTPUT_DIR)/*bowtie_final.aln

######################################
## Hi-C map processing 
##
######################################

## Only keep alignments which have both ends uniquely mapped to the genome
## Inputs : final.aln files (R1 + R2)
## Outputs : out files
separateAlignment_all:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Filter out reads alignement ..." >> $(LOGFILE)
	perl $(SCRIPTS)/separatePEalignment.pl -a $(BOWTIE2_FINAL_OUTPUT_DIR)/$(RES_FILE_NAME)_R1.final.aln -b $(BOWTIE2_FINAL_OUTPUT_DIR)/$(RES_FILE_NAME)_R2.final.aln -g $(ORGANISM) -o $(BOWTIE2_FINAL_OUTPUT_DIR)/$(RES_FILE_NAME)

## Assign alignments to regions segmented by HindIII sites
Mapped2HiCFragments_all:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Assign alignments to HindIII sites ..." >> $(LOGFILE)
	perl $(SCRIPTS)/overlapMapped2HiCFragments.pl -f1 $(GENOME_FRAGMENT) -f2 $(GENOME_FRAGMENT)  -m1 $(BOWTIE2_FINAL_OUTPUT_DIR)/$(RES_FILE_NAME).$(ORGANISM).1.out -m2 $(BOWTIE2_FINAL_OUTPUT_DIR)/$(RES_FILE_NAME).$(ORGANISM).2.out > $(LOGS_DIR)/overlapRS.log
	mv -f $(RES_FILE_NAME).$(ORGANISM)* $(DATA_DIR)
	gunzip -f $(DATA_DIR)/$(RES_FILE_NAME).$(ORGANISM)*.gz

## assign segmented regions to defined bins
AssignRead2bin_all:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Generate binned matrix files ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)), $(SCRIPTS)/assignRead2bins.sh -c $(CURDIR)/$(CONFIG_FILE) -i $(DATA_DIR)/$(RES_FILE_NAME).$(ORGANISM).interaction -g $(ANNOT_DIR)/chrom.sizes -b $(BSIZE);)

plot_FragmentInfo:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Plot Fragment Informations ..." >> $(LOGFILE)
	$(R_PATH)/R --no-save CMD BATCH "--args picDir='$(PIC_DIR)' dataDir='$(DATA_DIR)' sampleName='$(RES_FILE_NAME)'" $(SCRIPTS)/plotHiCFragment.R $(LOGS_DIR)/plotHiCFragment.Rout

matrix2RData:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Generate RData files ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)), mkdir -p $(MATRIX_DIR)/$(BSIZE);$(R_PATH)/R --no-save CMD BATCH "--args matDir='$(MATRIX_DIR)/$(BSIZE)' obj='$(RES_FILE_NAME_OBJ)_$(BSIZE)' cpu='$(N_CPU)' rdataDir='$(RDATA_DIR)' org='$(ORGANISM)'" $(SCRIPTS)/matrix2RData.R $(LOGS_DIR)/matrix2RData.Rout;)

######################################
## Normalization
##
######################################

ICEnorm:
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Run ICE Normalization ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),$(R_PATH)/R --no-save CMD BATCH "--args rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE).RData' norm='ICE' cpu='$(N_CPU)' outDir='$(RDATA_DIR)'" $(SCRIPTS)/hicNorm.R $(LOGS_DIR)/hicNorm.Rout;)
	@echo "Generate Matrix files from RData ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),mkdir -p $(ICED_MATRIX_DIR)/$(BSIZE); $(R_PATH)/R --no-save CMD BATCH "--args matDir='$(ICED_MATRIX_DIR)/$(BSIZE)' cpu='$(N_CPU)' rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_iced.RData' org='$(ORGANISM)'" $(SCRIPTS)/rData2matrix.R $(LOGS_DIR)/rData2matrix.Rout;)


LGFnorm: 
	@echo "--------------------------------------------" >> $(LOGFILE)
	@date >> $(LOGFILE)
	@echo "Run LGF Normalization ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),$(R_PATH)/R --no-save CMD BATCH "--args rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_annot.RData' norm='LGF' cpu='$(N_CPU)' outDir='$(RDATA_DIR)'" $(SCRIPTS)/hicNorm.R $(LOGS_DIR)/hicNorm.Rout;)
	@echo "Generate Matrix files from RData ..." >> $(LOGFILE)
	$(foreach BSIZE,$(subst $(comma),$(space),$(BIN_SIZE)),mkdir -p $(LGF_MATRIX_DIR)/$(BSIZE); $(R_PATH)/R --no-save CMD BATCH "--args matDir='$(LGF_MATRIX_DIR)/$(BSIZE)' cpu='$(N_CPU)' rdata='$(RDATA_DIR)/$(RES_FILE_NAME_OBJ)_$(BSIZE)_annot_lgf.RData' org='$(ORGANISM)'" $(SCRIPTS)/rData2matrix.R $(LOGS_DIR)/rData2matrix.Rout;)

######################################
## Clean
##
######################################

clean:
	/bin/rm -f $(BOWTIE2_GLOBAL_OUTPUT_DIR)/*.sam
	/bin/rm -f $(BOWTIE2_LOCAL_OUTPUT_DIR)/*.sam

reset:
	/bin/rm -rf bowtie_results hic_results $(LOGS_DIR)/*
