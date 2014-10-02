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

opts="-v"
if [[ "${ALL_OUTPUT}" -eq "1" ]]; then opts=$opts" -a"; fi
if [[ "${MIN_INSERT_SIZE}" -ge "0" ]]; then opts=$opts" -s ${MIN_INSERT_SIZE}"; fi
if [[ "${MAX_INSERT_SIZE}" -ge "0" ]]; then opts=$opts" -l ${MAX_INSERT_SIZE}"; fi

for r in $(get_files_for_overlap 1)
do
    sample_dir=$(get_sample_dir ${r})
    datadir=${MAPC_OUTPUT}/data/${sample_dir}
    mkdir -p ${datadir}
    
    cmd="python ${SCRIPTS}/overlapMapped2HiCFragments.py ${opts} -f ${GENOME_FRAGMENT} -r ${r} -o ${datadir}"
    exec_cmd $cmd
done
