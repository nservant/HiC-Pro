#!/bin/bash

## HiC-Pro
## Copyleft 2015-2018 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Combine two steps alignment files
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
## Combine both Bowtie2 mapping
##
mapping_combine()
{
    local sample_dir="$1"
    local file="$2"
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.bam//')
    local tmp_prefix=$(basename $prefix)

    mkdir -p ${BOWTIE2_FINAL_OUTPUT_DIR}/${sample_dir}    

    ## Set a default for legacy config files that do not have SORT_RAM set
    if [[ "${SORT_RAM}" == "" ]]; then
       SORT_RAM="768"
    fi

    ## Divide the SORT_RAM by the number of Cpus
    SORT_RAM=$(echo ${SORT_RAM} | sed -e 's/M$//')
    SORT_RAM=$(( ${SORT_RAM}/${N_CPU} ))
    
    
    ## Merge local and global alignment
    if [[ -e ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.bam && -e ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam ]]; then
	
	cmd="${SAMTOOLS_PATH}/samtools merge -@ ${N_CPU} -n -f ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam ${BOWTIE2_GLOBAL_OUTPUT_DIR}/${prefix}.bwt2glob.bam ${BOWTIE2_LOCAL_OUTPUT_DIR}/${prefix}.bwt2glob.unmap_bwt2loc.bam"
	exec_cmd $cmd 2>&1

        ## Sort merge file. In theory, should be perform by "merge -n", but do not work in some cases ... depending on read name ?
	cmd="${SAMTOOLS_PATH}/samtools sort -@ ${N_CPU} -m ${SORT_RAM}M -n -T ${TMP_DIR}/$tmp_prefix -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.sorted.bam ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
        exec_cmd $cmd 2>&1
    
	cmd="mv ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.sorted.bam ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam"
        exec_cmd $cmd 2>&1
    
    elif [[ ${LIGATTION_SITE} == "" ]]; then
	
	cmd="ln -f -s ../../bwt2_global/${prefix}.bwt2glob.bam ${BOWTIE2_FINAL_OUTPUT_DIR}/${prefix}.bwt2merged.bam "
	exec_cmd $cmd 2>&1
    else
	die "Error - Mapping files not found"
    fi
}

## Combine local and global alignments
for r in $(get_sam_for_combine)
do
    R1=$r
    R2=$(echo $r | get_R2)
    sample_dir=$(get_sample_dir $r)

    ## Logs
    ldir=${LOGS_DIR}/${sample_dir}
    echo "Logs: $ldir/mapping_combine.log"
    
    mapping_combine $sample_dir $R1 >> ${ldir}/mapping_combine.log &
    mapping_combine $sample_dir $R2 >> ${ldir}/mapping_combine.log &
    
    wait
done
