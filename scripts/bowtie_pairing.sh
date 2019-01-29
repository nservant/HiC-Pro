#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Launcher of mergeSAM script
## Merge R1 and R2 SAM files into one paired-end SAM file
##

dir=$(dirname $0)

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

CONF=$conf_file . $dir/hic.inc.sh

##
## Merge SE file into one PE file and filter the reads
##
merge_pairs()
{
    local sample_dir="$1"
    local file_r1="$2"
    local file_r2="$3"
    local ldir="$4"
    
    local prefix_r1=$(echo ${sample_dir}/$(basename $file_r1) | sed -e 's/.bwt2merged.bam//')
    local prefix_r2=$(echo ${sample_dir}/$(basename $file_r2) | sed -e 's/.bwt2merged.bam//')
    local prefix_out=$(echo $prefix_r1 | get_pairs)

    ## Create output folder if not available
    mkdir -p ${BOWTIE2_FINAL_OUTPUT_DIR}/${sample_dir}
    
    ## Merge two SAM files into 1 paired-end SAM file / removed unmapped and multihits reads
    OPTS="-q ${MIN_MAPQ} -t -v"
    if [[ ${RM_SINGLETON} == 0 ]]; then
	OPTS=$OPTS" -s"
    fi
    if [[ ${RM_MULTI} == 0 ]]; then
	OPTS=$OPTS" -m"
    fi

    #cmd="${PYTHON_PATH}/python ${SCRIPTS}/mergeSAM.py ${OPTS} -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r1}.bwt2merged.bam -r ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r2}.bwt2merged.bam -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam > ${ldir}/mergeSAM.log"
    cmd="${PYTHON_PATH}/python ${SCRIPTS}/mergeSAM.py ${OPTS} -f ${file_r1} -r ${file_r2} -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam"
    exec_cmd $cmd 
}

##
## Tag reads according to their SNP information
##
tag_allele_spe()
{
    local bam_paired="$1"
    local vcf_file="$2"
    local asout=$(echo ${sample_dir}/$(basename $r) | sed -e 's/.bwt2pairs.bam/.bwt2pairs_allspe.bam/')

    if [ -e ${vcf_file} ]; then
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/markAllelicStatus.py -s ${vcf_file} -v -r -i ${bam_paired} -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${asout}"
	exec_cmd $cmd
	cmd="mv ${BOWTIE2_FINAL_OUTPUT_DIR}/${asout} $bam_paired"
	exec_cmd $cmd
    else
	die "Error - VCF file not found"
    fi
}

## Combine R1/R2 tags in a single BAM file
for r in $(get_sam_for_merge)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    ## Logs
    ldir=${LOGS_DIR}/${sample_dir}
    mkdir -p ${ldir}
    echo "Logs: ${ldir}/mergeSAM.log"
    merge_pairs $sample_dir $R1 $R2 &> ${ldir}/mergeSAM.log
done

## Add allele specific tag if specified
if [[ ${ALLELE_SPECIFIC_SNP} != "" ]]; then
    as_file=`abspath $ALLELE_SPECIFIC_SNP`
    if [[ ${as_file} == "" || ! -f ${as_file} ]]; then
	as_file=$ANNOT_DIR/${ALLELE_SPECIFIC_SNP}
	if [[ ! -f $as_file ]]; then
	    echo "ALLELE_SPECIFIC_SNP not found. Cannot process alignment file"
	    exit 1
	fi
    fi

    for r in $(get_paired_bam)
    do
	sample_dir=$(get_sample_dir $r)
	ldir=${LOGS_DIR}/${sample_dir}
	echo "Logs: ${ldir}/markAllelicStatus.log"
	tag_allele_spe $r ${as_file} &> ${ldir}/markAllelicStatus.log
    done
fi
