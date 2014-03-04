#!/bin/bash
## Nicolas Servant
## Eric Viara updated 2014-04-28
##

dir=$(dirname $0)

. $dir/hic.inc.sh

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

read_config $ncrna_conf

mapping_combine()
{
    local sample_dir="$1"
    local file="$2"
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.aln//')

    echo ${prefix} >> ${LOGFILE}

    mkdir -p ${BOWTIE2_FINAL_OUTPUT_DIR}/${sample_dir}

    cmd="${SCRIPTS}/mergeBwt2GlobLoc.pl -a ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.aln -b ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.aln -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}_bwt2merged.aln"
    exec_cmd $cmd
}

for r in $(get_aln_for_merge)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    mapping_combine $sample_dir $R1 &
    mapping_combine $sample_dir $R2 &

    wait

    if [ 0 = 1 ]; then
    echo ${prefix} >> ${LOGFILE}

    perl ${SCRIPTS}/mergeBwt2GlobLoc.pl -a ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.aln -b ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.aln -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}_bwt2merged.aln;
    fi
done
