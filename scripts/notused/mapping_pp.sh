#!/bin/bash
## Nicolas Servant
## Institut Curie

## Mapping post_process

dir=$(dirname $0)

. $dir/hic.inc.sh

## Usage
function usage {
    echo -e "Usage : ./mapping_pp.sh"
    echo -e "-u"" Generate unmapped fastq file"
    echo -e "-l"" local"
    echo -e "-h"" <help>"
    exit
}

unmap=0
local=0

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-u) unmap=1; shift;;
	(-l) local=1; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

read_config $ncrna_conf

mapping_pp()
{
    local file="$1"
    local output_aln="$2"
    local cmd="perl ${SCRIPTS}/BWT2output2Novoformat.pl -a ${file} -o ${output_aln}"
    exec_cmd "$cmd"
    
    ## Get unmapped reads
    if [[ $unmap == 1 ]]; then
	local output_unmapped=`echo ${file} | sed -e 's/bam/unmap.fastq/'`
	
	cmd="perl ${SCRIPTS}/extractUndefineSeq.pl -a $output_aln -o ${output_unmapped}"
	exec_cmd "$cmd"
    fi
}

echo "mapping_fastq_pp local $local"
for r in $(get_bam_for_pp $local)
do
    R1=$r
    R2=$(echo $r | get_R2)
    output_aln1=$(echo $R1 | sed -e 's/bam/aln/')
    output_aln2=$(echo $R2 | sed -e 's/bam/aln/')
    
    mapping_pp $R1 $output_aln1 &
    mapping_pp $R2 $output_aln2 &

    wait

    if [ 0 = 1 ]; then
    ## Post-processing
    cmd="perl ${SCRIPTS}/BWT2output2Novoformat.pl -a ${r} -o ${output_aln}"
    echo $cmd
    eval $cmd
    
    ## Get unmapped reads
    if [[ $unmap == 1 ]]; then
    	output_unmapped=`echo ${r} | sed -e 's/bam/unmap.fastq/'`
	
    	cmd="perl ${SCRIPTS}/extractUndefineSeq.pl -a $output_aln -o ${output_unmapped}"
    	echo $cmd
    	eval $cmd
     fi
    fi
done
