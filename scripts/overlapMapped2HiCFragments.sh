#!/bin/bash
## Eric Viara 2014-04-30
##

dir=$(dirname $0)

. $dir/hic.inc.sh

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

read_config $ncrna_conf

for r in $(get_files_for_overlap 1)
do
    R1=$r
    R2=$(echo $R1 | sed -e 's/1\.out/2.out/')
    file2move=$(basename $(echo $R1 | sed -e 's/\.1\.out//'))
    sample_dir=$(get_sample_dir $R1)

    cmd="perl ${SCRIPTS}/overlapMapped2HiCFragments.pl -f1 ${GENOME_FRAGMENT} -f2 ${GENOME_FRAGMENT}  -m1 $R1 -m2 $R2"
    exec_cmd $cmd

    datadir=${MAPC_OUTPUT}/data/${sample_dir}
    mkdir -p ${datadir}

    gunzip -f ${file2move}*.gz
    mv ${file2move}* ${datadir}
done
