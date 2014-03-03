#!/bin/bash
## Nicolas Servant
## Institut Curie

## Mapping post_process

## Usage
function usage {
    echo -e "Usage : ./mapping_pp.sh"
    echo -e "-i"" <input directory>"
    echo -e "-u"" Generate unmapped fastq file"
    echo -e "-h"" <help>"
    exit
}


################### Initialize ###################
unmap=0
set -- $(getopt i:u "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-i) input_dir=$2; shift;;
	(-u) unmap=1; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

if [ ! -d $input_dir ]
then
    echo "$input_dir" is not a directory
    usage
    exit
fi

## Alignment Post-process
SCRIPTS=`dirname $0`

for r in ${input_dir}/*.bam
do
    output_aln=`echo ${r} | sed -e 's/bam/aln/'`
    
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
done


