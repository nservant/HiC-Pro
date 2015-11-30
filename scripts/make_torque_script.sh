#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Create PBS Torque files
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
    inputfile=inputfiles_${PBS_SUFFIX}.txt
    get_hic_files $RAW_DIR .fastq | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
elif [[ $MAKE_OPTS == *"proc_hic"* ]]
then
    inputfile=inputfiles_${PBS_SUFFIX}.txt
    get_hic_files $RAW_DIR .bam | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
fi

## Paralelle Implementation
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"mapping"* || $MAKE_OPTS == *"proc_hic"* ]]
then
    make_target="all_qsub"
    ## Remove per sample steps
    if [[ $MAKE_OPTS != "" ]]; then 
	make_target=$(echo $MAKE_OPTS | sed -e 's/,/ /g'); 
	make_target=$(echo $make_target | sed -e 's/merge_persample//g');
	make_target=$(echo $make_target | sed -e 's/build_contact_maps//g');
	make_target=$(echo $make_target | sed -e 's/ice_norm//g');
        make_target=$(echo $make_target | sed -e 's/quality_checks//g');
    fi
 
    ## step 1 - parallel
    torque_script=HiCPro_step1_${PBS_SUFFIX}.sh
    PPN=$(( ${N_CPU} * 2))
    cat > ${torque_script} <<EOF
#!/bin/bash
#PBS -l nodes=1:ppn=${PPN},mem=${PBS_MEM},walltime=${PBS_WALLTIME}
#PBS -M ${PBS_MAIL}
#PBS -m ae
#PBS -j eo
#PBS -N HiCpro_s1_${PBS_SUFFIX}
#PBS -q ${PBS_QUEUE}
#PBS -V
#PBS -t 1-$count

cd \$PBS_O_WORKDIR

FASTQFILE=\$PBS_O_WORKDIR/$inputfile; export FASTQFILE
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${torque_script}

    ## User message
    echo "The following command will launch the parallel workflow through $count torque jobs:"
    echo qsub ${torque_script}
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

    torque_script_s2=HiCPro_step2_${PBS_SUFFIX}.sh
    cat > ${torque_script_s2} <<EOF
#!/bin/bash
#PBS -l nodes=1:ppn=1,mem=${PBS_MEM},walltime=${PBS_WALLTIME}
#PBS -M ${PBS_MAIL}
#PBS -m ae
#PBS -j eo
#PBS -N HiCpro_s2_${PBS_SUFFIX}
#PBS -q ${PBS_QUEUE}
#PBS -V

cd \$PBS_O_WORKDIR
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${torque_script_s2}

    ## User message
    echo "The following command will merge the processed data and run the remaining steps per sample:"
    echo qsub ${torque_script_s2}
fi

