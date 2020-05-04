#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Guipeng Li, Nicolas Servant 
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## Adapted for LSF (specifically, the PMACS HPC Cluster); asrini@upenn.edu 2016-02-29

##
## Create LSF files
##

dir=$(dirname $0)

usage()
{
    echo "usage: $0 -c CONFIG [-s STEP]"
}

MAKE_OPTS=""

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(-s) MAKE_OPTS=$2; shift;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  suffix=$1; break;;
    esac
    shift
done

if [ -z "$conf_file" ]; then usage; exit 1; fi

CONF=$conf_file . $dir/hic.inc.sh
unset FASTQFILE

## Define input files
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"mapping"* ]]
then
    inputfile=inputfiles_${JOB_NAME}.txt
    ifq=$(get_hic_files $RAW_DIR .fq)
    ifastq=$(get_hic_files $RAW_DIR .fastq)
    echo -e "$ifq\n$ifastq" | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
elif [[ $MAKE_OPTS == *"proc_hic"* ]]
then
    inputfile=inputfiles_${JOB_NAME}.txt
    get_hic_files $RAW_DIR .bam | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
fi

## Paralelle Implementation
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"mapping"* || $MAKE_OPTS == *"proc_hic"* ]]
then
    make_target="all_sub"
    ## Remove per sample steps
    if [[ $MAKE_OPTS != "" ]]; then 
	make_target=$(echo $MAKE_OPTS | sed -e 's/,/ /g'); 
	make_target=$(echo $make_target | sed -e 's/merge_persample//g');
	make_target=$(echo $make_target | sed -e 's/build_contact_maps//g');
	make_target=$(echo $make_target | sed -e 's/ice_norm//g');
    make_target=$(echo $make_target | sed -e 's/quality_checks//g');
    fi
 
    ## step 1 - parallel
    lsf_script=HiCPro_step1_${JOB_NAME}.sh
    cat > ${lsf_script} <<EOF
#!/bin/bash
#BSUB -M ${JOB_MEM}
#BSUB -W ${JOB_WALLTIME}
#BSUB -N 
#BSUB -u ${JOB_MAIL}
#BSUB -J HiCpro_s1_${JOB_NAME}[1-$count]
#BSUB -e HiCpro_s1_${JOB_NAME}.%J.e
#BSUB -o HiCpro_s1_${JOB_NAME}.%J.o
#BSUB -q ${JOB_QUEUE}
#BSUB -n ${N_CPU}
EOF

    if [[ $count -gt 1 ]]; then
	echo -e "#BSUB -J HiCpro_s1_${JOB_NAME}[1-$count]" >> ${lsf_script} 
    fi

    cat >> ${lsf_script} <<EOF
FASTQFILE=$inputfile; export FASTQFILE
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${lsf_script}

    ## User message
    echo "The following command will launch the parallel workflow through $count lsf jobs:"
    echo "bsub < ${lsf_script}"
fi    


## Per sample Implementation
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"build_contact_maps"* || $MAKE_OPTS == *"ice_norm"* || $MAKE_OPTS == *"quality_checks"* ]]
then
    make_target="all_persample"
    ## Remove parallele mode
    if [[ $MAKE_OPTS != "" ]]; 
    then 
	make_target=$(echo $MAKE_OPTS | sed -e 's/,/ /g'); 
	make_target=$(echo $make_target | sed -e 's/mapping//g');
	make_target=$(echo $make_target | sed -e 's/proc_hic//g');
    fi

    lsf_script_s2=HiCPro_step2_${JOB_NAME}.sh
    cat > ${lsf_script_s2} <<EOF
#!/bin/bash
#BSUB -M ${JOB_MEM}
#BSUB -W ${JOB_WALLTIME}
#BSUB -N 
#BSUB -u ${JOB_MAIL}
#BSUB -J HiCpro_s2_${JOB_SUFFIX}
#BSUB -e HiCpro_s2_${JOB_SUFFIX}.%J.e
#BSUB -o HiCpro_s2_${JOB_SUFFIX}.%J.o
#BSUB -q ${JOB_QUEUE}
#BSUB -n ${N_CPU}

make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${lsf_script_s2}

    ## User message
    echo "The following command will merge the processed data and run the remaining steps per sample:"
    echo "bsub < ${lsf_script_s2}"
fi

