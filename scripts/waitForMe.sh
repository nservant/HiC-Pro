#!/bin/bash


################### Initialize ###################
mode='global'
samplelist=""
set -- $(getopt c:p:l "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-l) mode='local'; shift;;
	(-p) samplelist=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done


sampe_jobname=`echo "myscript.sh'" | qsub -m abe $PBS_OUTPUT -e $PBS_ERROR -N ${SAMPLENAME}_SAMPE  -q batch -l nodes=1:ppn=8,mem=10gb`

sampe_jobid=${sampe_jobname%%.*}

JOB2WAIT="$PBS_OUTPUT/${SAMPLENAME}_SAMPE.o$sampe_jobid"

if [ -e ${JOB2WAIT} ];then 

	 STATUS="done"
else STATUS="running"
fi

while [ $STATUS != "done" ]
do
	sleep 2
	if [ -e ${JOB2WAIT} ];then

	     STATUS="done"
	else STATUS="running"
	fi

done

