#!/bin/bash
## Nicolas Servant updated 2014-08-13
## Eric Viara updated 2014-04-28
##

dir=$(dirname $0)
##. $dir/hic.inc.sh

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

##read_config $CONF
CONF=$conf_file . $dir/hic.inc.sh

##
## Merge SE file into one PE file and filter the reads
##
merge_pairs()
{
    local sample_dir="$1"
    local file_r1="$2"
    local file_r2="$3"

    local prefix_r1=$(echo ${sample_dir}/$(basename $file_r1) | sed -e 's/.bwt2merged.sam//')
    local prefix_r2=$(echo ${sample_dir}/$(basename $file_r2) | sed -e 's/.bwt2merged.sam//')
    local prefix_out=$(echo $prefix_r1 | get_pairs)

    ## Merge two SAM files into 1 paired-end SAM file / removed unmapped and multihits reads
    cmd="${SCRIPTS}/mergeSAM.pl -u -m -q ${MIN_MAPQ} -c ${CUT_SITE_5OVER} -v -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r1}.bwt2merged.sam -r ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_r2}.bwt2merged.sam -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.sam > ${LOGS_DIR}/"$(basename ${prefix_out})"_merge.log"
    exec_cmd $cmd

    ## Generate BAM file
    cmd="${SAMTOOLS_PATH}/samtools view -bS ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.sam > ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix_out}.bwt2pairs.bam"
    exec_cmd $cmd
}

## Combine R1/R2 tags in a single BAM file
for r in $(get_sam_for_merge)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    merge_pairs $sample_dir $R1 $R2 &
done
