#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Mapping statistics for Hi-C data
##

dir=$(dirname $0)

## Usage
function usage {
    echo -e "Usage : ./mapping_stat.sh"
    echo -e "-i"" <input directory>"
    echo -e "-c"" <config>"
    echo -e "-o"" <output directory/prefix>"
    echo -e "-h"" <help>"
    exit
}

mode=global
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-i) input_dir=$2; shift;;
	(-o) output_dir=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

CONF=$conf_file . $dir/hic.inc.sh

##
## Mapping stat
##
mapping_stat(){

    local sample_dir="$1"
    local file="$2"
    local tag="$3"
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.bam//')

    cmd="${SAMTOOLS_PATH}/samtools view -c ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
    echo $cmd >&2
    tot_reads=`exec_ret $cmd`
    cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
    echo $cmd >&2
    map_reads=`exec_ret $cmd`
    cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.bam"
    echo $cmd >&2
    gmap_reads=`exec_ret $cmd`
    if [[ -e ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam ]]; then
	cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam"
	echo $cmd >&2
	lmap_reads=`exec_ret $cmd`
    else
	lmap_reads=0
    fi
    echo "## HiC-Pro Mapping Statistics"
    echo "## $prefix.mapstat"
    echo -e "total_${tag}\t$tot_reads"
    echo -e "mapped_${tag}\t$map_reads"
    echo -e "global_${tag}\t$gmap_reads"
    echo -e "local_${tag}\t$lmap_reads"
}

for r in $(get_global_aln_for_stats ${mode})
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir ${r})

    ## logs
    ldir=${LOGS_DIR}/${sample_dir}
    
    R_STAT1=$(get_stat_file $R1)
    R_STAT2=$(get_stat_file $R2)

    echo "Logs: ${ldir}/mapping_stats.log"
    mapping_stat $sample_dir $R1 R1 1> $R_STAT1 2>> ${ldir}/mapping_stats.log &
    mapping_stat $sample_dir $R2 R2 1> $R_STAT2 2>> ${ldir}/mapping_stats.log &

    wait
done
