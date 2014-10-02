#!/bin/bash
## Nicolas Servant
##

OUTPUT="sample_list.csv"

################### Initialize ###################
set -- $(getopt i: "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-i) RAW_DIR=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

if [ -e $OUTPUT ];then
    rm $OUTPUT
fi

NB_SAMPLE=0
for r in ${RAW_DIR}/*.fastq
do
    ((NB_SAMPLE++))
    echo $NB_SAMPLE";"${r} >> $OUTPUT
done

echo $NB_SAMPLE