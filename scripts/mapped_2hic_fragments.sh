#!/bin/bash
## Eric Viara 2014-04-30
##

dir=$(dirname $0)

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

##read_config $ncrna_conf
CONF=$conf_file . $dir/hic.inc.sh

opts="-v"
if [[ "${GET_ALL_INTERACTION_CLASSES}" -eq "1" ]]; then opts=$opts" -a"; fi
if [[ "${GET_PROCESS_SAM}" -eq "1" ]]; then opts=$opts" -S"; fi
if [[ "${MIN_INSERT_SIZE}" -ge "0" && "${MIN_INSERT_SIZE}" -ne "" ]]; then opts=$opts" -s ${MIN_INSERT_SIZE}"; fi
if [[ "${MAX_INSERT_SIZE}" -ge "0" && "${MIN_INSERT_SIZE}" -ne "" ]]; then opts=$opts" -l ${MAX_INSERT_SIZE}"; fi

GENOME_FRAGMENT_FILE=`abspath $GENOME_FRAGMENT`
if [[ ! -e $GENOME_FRAGMENT_FILE ]]; then
    GENOME_FRAGMENT_FILE=$ANNOT_DIR/$GENOME_FRAGMENT
    if [[ ! -e $GENOME_FRAGMENT_FILE ]]; then
	echo "$GENOME_FRAGMENT not found. Exit"
	exit -1
    fi
fi

for r in $(get_files_for_overlap)
do
    sample_dir=$(get_sample_dir ${r})
    datadir=${MAPC_OUTPUT}/data/${sample_dir}
    mkdir -p ${datadir}
    
    cmd="python ${SCRIPTS}/mapped_2hic_fragments.py ${opts} -f ${GENOME_FRAGMENT_FILE} -r ${r} -o ${datadir}"
    exec_cmd $cmd

    ## Valid pairs are already sorted
    outVALID=`basename ${r} | sed -e 's/.sam$/.validPairs/'`
    outSAM=`basename ${r} | sed -e 's/.sam$/_interaction.sam/'`
    sortBAM=`basename ${r} | sed -e 's/.sam$/_interaction/'`
    
    echo "## Sorting valid interaction file ..."
    sort -k2,2V -k3,3n -k5,5V -k6,6n -T ${TMP_DIR} -o ${datadir}/${outVALID} ${datadir}/${outVALID} 

    echo "## Creating BAM and index files ..."
    if [ -f ${datadir}/${outSAM} ]
    then 
	cmd="${SAMTOOLS_PATH}/samtools view -bS ${datadir}/${outSAM} | ${SAMTOOLS_PATH}/samtools sort - ${datadir}/${sortBAM}"
	exec_cmd $cmd
	cmd="${SAMTOOLS_PATH}/samtools index ${datadir}/${sortBAM}.bam"
	exec_cmd $cmd
    fi
done
echo "## done !"
