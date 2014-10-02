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

echo "separatePE"

for r in $(get_aln_for_separatePE)
do
    R1=$r
    R2=$(echo $r | get_R2)
    outfile=$(basename $(echo $R1 | sed -e 's/_R1.*$//'))
    sample_dir=$(get_sample_dir $R1)

    cmd="perl ${SCRIPTS}/separatePEalignment.pl -a $R1 -b $R2 -g ${ORGANISM} -o ${BOWTIE2_FINAL_OUTPUT_DIR}/${sample_dir}/${outfile}"
    exec_cmd $cmd
done
