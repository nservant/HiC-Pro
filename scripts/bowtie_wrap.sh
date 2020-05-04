#!/bin/bash

## HiC-Pro
## Copyright (c) 2015-2018 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## Bowtie2 wrapper
## Step1 mapping
end_to_end_align()
{
    local odir="$1"
    local infile="$2"
    local unmap="$3"
    local filtunmap="$4"
    local ldir="$5"
    
    ## Check mapping options
    if [[ -z ${BOWTIE2_GLOBAL_OPTIONS} ]]; then
	echo "Error: Mapping step1 options not defined. Exit"
	exit -1
    fi
    
    ## Output
    prefix=$(basename ${infile} | sed -e 's/.fastq\(.gz\)*$//' -e 's/.fq\(.gz\)*$//')
    
    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_GLOBAL_OPTIONS=${BOWTIE2_GLOBAL_OPTIONS}" --un ${odir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
    fi
    
    ## Check for gz files
    if [[ ! -e $infile && -e $infile.gz ]]; then
	infile=$infile.gz
    fi
    
    ## Run bowtie
    if [[ $N_CPU -lt 2 ]]; then
	echo -e "Warning : HiC-Pro need at least 2 CPUs to run the mapping !!"
	bwt_cpu=1
    else
	bwt_cpu=$(( $N_CPU / 2 ))
    fi
    
    echo "##HiC-Pro mapping" > ${ldir}/${prefix}_bowtie2.log
    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_GLOBAL_OPTIONS} --rg-id BMG --rg SM:${prefix} -p ${bwt_cpu} -x ${BOWTIE2_IDX} -U ${infile} 2>> ${ldir}/${prefix}_bowtie2.log"
    if [[ $filtunmap == 1 ]]; then
	cmd=$cmd"| ${SAMTOOLS_PATH}/samtools view -F 4 -bS - > ${odir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam"
    else
	cmd=$cmd"> ${odir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.bam"
    fi
   
    exec_cmd $cmd 
}

## Bowtie2
## Trimmed reads Alignment
cut_and_align()
{
    local odir="$1"
    local infile="$2"
    local unmap="$3"
    local ldir="$4"
    
    ## Check mapping options
    if [[ -z ${BOWTIE2_LOCAL_OPTIONS} ]]; then
	echo "Mapping step2 options not defined. Exit"
	exit -1
    fi
    ## Output
    prefix=$(basename ${infile} | sed -e 's/.fastq\(.gz\)*$//' -e 's/.fq\(.gz\)*$//')

    ## Starts trimming reads from the ligation site
    tfile=`basename $infile | sed -e s/.fastq$/_trimmed.fastq/`
    ##cmd="${SCRIPTS}/cutsite_trimming --fastq $infile --cutsite ${LIGATION_SITE} --out ${odir}/$tfile --rmuntrim > ${ldir}/${prefix}_readsTrimming.log 2>&1"
    cmd="${SCRIPTS}/cutsite_trimming --fastq $infile --cutsite ${LIGATION_SITE} --out ${odir}/$tfile  > ${ldir}/${prefix}_readsTrimming.log 2>&1"
    exec_cmd "$cmd" 
    
    ## Unmapped reads
    if [[ $unmap == 1 ]]; then
	BOWTIE2_LOCAL_OPTIONS=${BOWTIE2_LOCAL_OPTIONS}" --un ${odir}/${prefix}_${REFERENCE_GENOME}.bwt2glob.unmap.fastq"
    fi

    if [[ $N_CPU -lt 2 ]]; then
        echo -e "Warning : HiC-Pro need at least 2 CPUs to run the mapping !!"
        bwt_cpu=1
    else
        bwt_cpu=$(( $N_CPU / 2 ))
    fi
       
    ## Run bowtie
    echo "##HiC-Pro mapping" > ${ldir}/${prefix}_bowtie2.log

    cmd="${BOWTIE2_PATH}/bowtie2 ${BOWTIE2_LOCAL_OPTIONS} --rg-id BML --rg SM:${prefix} -p ${bwt_cpu} -x ${BOWTIE2_IDX} -U ${odir}/${tfile} 2>> ${ldir}/${prefix}_bowtie2.log | ${SAMTOOLS_PATH}/samtools view -bS - > ${odir}/${prefix}_bwt2loc.bam"
    exec_cmd "$cmd" 
}


## Init
dir=$(dirname $0)
MODE='step1'
UNMAP=0

## Get args
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-m) MODE=$2; shift;;
	(-p) FASTQLIST=$2; shift;;
	(-u) UNMAP=1; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done


## Read configuration files
CONF=$conf_file . $dir/hic.inc.sh

## If the step 2 is not run - do not filter the BAM
if [[ -z $LIGATION_SITE ]]; then
    FILT_UNMAP=0
else
    FILT_UNMAP=1
fi

if [[ ${MODE} == 'step1' ]]; then
    for r in $(get_fastq_for_bowtie_global)
    do
	wasrun=1
       	R1=$r
	R2=$(echo $r | get_R2)
	sample_dir=$(get_sample_dir $r)

	## Ouput
	odir=${BOWTIE2_GLOBAL_OUTPUT_DIR}/${sample_dir}
	mkdir -p $odir
	
	## Logs
	ldir=${LOGS_DIR}/${sample_dir}
	mkdir -p ${ldir}

	echo "Logs: ${ldir}/mapping_step1.log"
	
	end_to_end_align ${odir} ${R1} ${UNMAP} ${FILT_UNMAP} ${ldir} >> ${ldir}/mapping_step1.log & 
	pid1=$!
	end_to_end_align ${odir} ${R2} ${UNMAP} ${FILT_UNMAP} ${ldir} >> ${ldir}/mapping_step1.log &
	pid2=$!
	
	wait $pid1 $pid2 || die "Error in reads alignment - Exit"
    done
    if [[ $wasrun != 1 ]]; then
	echo "Nothing to align ! Please check input files and R1/R2 extension." >&2
	exit 1
    fi
elif [[ ${MODE} == 'step2' ]]; then
    if [[ ! -z $LIGATION_SITE ]]; then
	for r in $(get_fastq_for_bowtie_local)
	do
	    wasrun=1
	    R1=$r
	    R2=$(echo $r | get_R2)	    	    
	    sample_dir=$(get_sample_dir $r)

	    ## Ouput
	    odir=${BOWTIE2_LOCAL_OUTPUT_DIR}/${sample_dir}
	    mkdir -p $odir

	    ## Logs
	    ldir=${LOGS_DIR}/${sample_dir}
	    mkdir -p ${ldir}

	    echo "Logs: ${ldir}/mapping_step2.log"
	    
	    cut_and_align ${odir} ${R1} ${UNMAP} ${ldir} >> ${ldir}/mapping_step2.log &
	    pid1=$!
	    cut_and_align ${odir} ${R2} ${UNMAP} ${ldir} >> ${ldir}/mapping_step2.log &
	    pid2=$!
	    
	    wait $pid1 $pid2 || die "Error in reads alignment - Exit"
	done
	if [[ $wasrun != 1 ]]; then
	    echo "Nothing to align ! Please check input files and R1/R2 extension." >&2
	    exit 1
	fi
    fi
else
    die "Error: Unknown mapping mode !"
fi
