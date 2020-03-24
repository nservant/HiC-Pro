#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Merge all valid interaction file and remove duplicates
## Split the interaction into G1/G2 files if specified
##
dir=$(dirname $0)

################### Initialize ###################

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

################### Read the config file ###################

##read_config $ncrna_conf
CONF=$conf_file . $dir/hic.inc.sh

################### Define Input Directory ###################
input_data_type=$(get_data_type)

if [[ $input_data_type == "valid" || $input_data_type == "allvalid" ]]
then
    IN_DIR=${RAW_DIR}
else
    IN_DIR=${MAPC_OUTPUT}/data/
fi

nbf=$(find -L ${IN_DIR} -mindepth 1 | wc -l)
if [[ $nbf == 0 ]]; then die "Error : empty ${IN_DIR} folder."; fi

################### Merge valid interaction files ###################
for RES_FILE_NAME in ${IN_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    
    ## out
    DATA_DIR=${MAPC_OUTPUT}/data/
    mkdir -p ${DATA_DIR}/${RES_FILE_NAME}
    STATS_DIR=${MAPC_OUTPUT}/stats/
    mkdir -p ${STATS_DIR}/${RES_FILE_NAME}


    ## Logs
    ldir=${LOGS_DIR}/${RES_FILE_NAME}
    mkdir -p ${ldir}
    logfile=${ldir}/merge_valid_interactions.log
    
    echo "Logs: ${logfile}"
    echo "## Merge valid interactions for ${RES_FILE_NAME} - ${input_data_type} - ${IN_DIR}..." > ${logfile}
  
    if [ -d ${IN_DIR}/${RES_FILE_NAME} ]; then

	## Check data dir is not empty
	nbf=$(find -L ${IN_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.validPairs" | wc -l)
	if [[ $nbf == 0 ]]; then die "Error : no valid interaction files found in ${DATA_DIR}/${RES_FILE_NAME}."; fi

	## Rm duplicates
	if [[ ${RM_DUP} == 0 || -z ${RM_DUP} ]]
	then
	    echo "## Do NOT remove duplicates ..." >> ${logfile}
	    cmd="cat ${IN_DIR}/${RES_FILE_NAME}/*.validPairs > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs"
	    exec_cmd $cmd >> ${logfile} 2>&1
	elif [[ ${RM_DUP} == 1 ]]; then
	     echo "## Remove duplicates ..." >> ${logfile}
	     cmd="LANG=en; sort -T ${TMP_DIR} -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -m ${IN_DIR}/${RES_FILE_NAME}/*.validPairs | \
awk -F\"\t\" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=\$2 || c2!=\$5 || s1!=\$3 || s2!=\$6){print;c1=\$2;c2=\$5;s1=\$3;s2=\$6}' > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs"
	    exec_cmd $cmd >> ${logfile} 2>&1
	fi

	allcount=$(cat  ${IN_DIR}/${RES_FILE_NAME}/*.validPairs | wc -l)
	allcount_rmdup=$(cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs | wc -l)
	nbdup=$(( allcount-allcount_rmdup ))

        ## merge stat file
	echo -e $RES_FILE_NAME"\t"$nbdup >> ${logfile}
        echo -e "valid_interaction\t"$allcount > ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat
        echo -e "valid_interaction_rmdup\t"$allcount_rmdup >> ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat
        awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} $2 == $5{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1}END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs >> ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat

       ## On Target                                                                                                                                                      
        if [[ ! -z ${CAPTURE_TARGET} ]]; then
            echo "## Select valid interactions from capture target ..." >> ${logfile}
	    if [[ ${REPORT_CAPTURE_REPORTER} == "0" ]]; then
		cap_opts="--cis";
	    fi
	    cmd="${PYTHON_PATH}/python ${SCRIPTS}/onTarget.py -i ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs -t ${CAPTURE_TARGET} $cap_opts -s ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat -v > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_ontarget.allValidPairs"
	    exec_cmd $cmd >> ${logfile} 2>&1
	fi
	
	## Allele specific analysis
	if [[ ! -z ${ALLELE_SPECIFIC_SNP} ]]; then
	    echo "## Split valid interactions for allele specific maps ..." >> ${logfile}
	    if [[ ! -z ${CAPTURE_TARGET} ]]; then
		cmd="${PYTHON_PATH}/python ${SCRIPTS}/split_valid_interactions.py -i ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_ontarget.allValidPairs -s ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_\
allValidPairs_assplit.stat -v"
	    else
		cmd="${PYTHON_PATH}/python ${SCRIPTS}/split_valid_interactions.py -i ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}.allValidPairs -s ${STATS_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_\
allValidPairs_assplit.stat -v"
	    fi
	    exec_cmd $cmd >> ${logfile} 2>&1
	fi
    fi
    wait
done

