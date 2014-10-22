#!/bin/bash
## Nicolas Servant
## Institut Curie

## Mapping statistics for Hi-C data

dir=$(dirname $0)

#. $dir/hic.inc.sh

## Usage
function usage {
    echo -e "Usage : ./mappingstat.sh"
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
	(-l) mode=local; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

##read_config $ncrna_conf
CONF=$conf_file . $dir/hic.inc.sh

##
## Mapping stat
##
mapping_stat(){

    local sample_dir="$1"
    local file="$2"
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.sam//')

    cmd="${SAMTOOLS_PATH}/samtools view -c ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
    tot_reads=`exec_ret $cmd`
    cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
    map_reads=`exec_ret $cmd`
    cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.bam"
    gmap_reads=`exec_ret $cmd`
    cmd="${SAMTOOLS_PATH}/samtools view -c -F 4 ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam"
    lmap_reads=`exec_ret $cmd`
    
    echo "## $prefix.mapstat"
    echo "$tot_reads total"
    echo "$map_reads mapped"
    echo "$gmap_reads global"
    echo "$lmap_reads local"
}

for r in $(get_aln_for_stats ${mode})
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir ${r})

    R_STAT1=$(get_stat_file $mode $R1)
    R_STAT2=$(get_stat_file $mode $R2)
    
    mapping_stat $sample_dir $R1 > $R_STAT1 &
    mapping_stat $sample_dir $R2 > $R_STAT2 &

    wait
done
