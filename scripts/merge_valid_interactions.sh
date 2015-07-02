#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Merge all valid interaction file and remove duplicates
##

dir=$(dirname $0)

##. $dir/hic.inc.sh

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

################### Define Variables ###################

DATA_DIR=${MAPC_OUTPUT}/data/
mkdir -p ${DATA_DIR}
nbf=$(find -L ${DATA_DIR} -mindepth 1 | wc -l)
if [[ $nbf == 0 ]]; then die "Error : empty ${DATA_DIR} folder."; fi

################### Combine Bowtie mapping ###################
for RES_FILE_NAME in ${DATA_DIR}/*
do
    RES_FILE_NAME=$(basename $RES_FILE_NAME)
    ## Logs
    LDIR=${LOGS_DIR}/${RES_FILE_NAME}
    mkdir -p ${LDIR}

    echo "## Merge valid interactions for ${RES_FILE_NAME}..." > ${LDIR}/merge_valid_interactions.log    
  
    if [ -d ${DATA_DIR}/${RES_FILE_NAME} ]; then

	## Check data dir is not empty
	nbf=$(find -L ${DATA_DIR}/${RES_FILE_NAME} -maxdepth 1 -name "*.validPairs" | wc -l)
	if [[ $nbf == 0 ]]; then die "Error : no valid interaction files found in ${DATA_DIR}/${RES_FILE_NAME}."; fi

	if [[ ${RM_DUP} == 0 ]]
	then
	    cat ${DATA_DIR}/${RES_FILE_NAME}/*.validPairs > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs
	else
	    echo "## Remove duplicates ..." >> ${LDIR}/merge_valid_interactions.log
	    allcount=$(cat  ${DATA_DIR}/${RES_FILE_NAME}/*.validPairs | wc -l)
	    sort -k2,2V -k3,3n -k5,5V -k6,6n -T ${TMP_DIR} -m ${DATA_DIR}/${RES_FILE_NAME}/*.validPairs | awk -F"\t" 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs
	    allcount_rmdup=$(cat ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs | wc -l)
	    nbdup=$(( allcount-allcount_rmdup ))
	    echo -e $RES_FILE_NAME"\t"$nbdup >> ${LDIR}/merge_valid_interactions.log

	    ## merge stat file
	    echo -e "valid_interaction\t"$allcount > ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat
	    echo -e "valid_interaction_rmdup\t"$allcount_rmdup >> ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat
	    awk '$2 == $5{cis=cis+1; if ($6-$3<=20000){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1}END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs >> ${DATA_DIR}/${RES_FILE_NAME}/${RES_FILE_NAME}_allValidPairs.mergestat
	fi
     fi
    wait
done

## make plots
## ${SCRIPTS}/make_plots.sh -c ${conf_file} -p "contacts" >> ${LOGFILE}
