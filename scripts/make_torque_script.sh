#!/bin/bash
## HiC-Pro
## Copyleft 2015 Institut Curie                               
## Author(s): Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##
## Create PBS Torque files
##

dir=$(dirname $0)

usage()
{
    echo "usage: $0 -c CONFIG"
}

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf_file=$2; shift;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  suffix=$1; break;;
    esac
    shift
done

if [ -z "$conf_file" ]; then usage; exit 1; fi

##read_config $conf_file

CONF=$conf_file . $dir/hic.inc.sh
unset FASTQFILE

fastqfile=fastqfile_${PBS_SUFFIX}.txt
get_hic_files $RAW_DIR .fastq | sed -e "s|$RAW_DIR||" -e "s|^/||" > $fastqfile
count=$(cat $fastqfile | wc -l)

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

FASTQFILE=\$PBS_O_WORKDIR/$fastqfile; export FASTQFILE
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt all_qsub 2>&1
EOF

chmod +x ${torque_script}

## step 2
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
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt build_contact_maps 2>&1
EOF

chmod +x ${torque_script_s2}

echo "Please run HiC-Pro in two steps :"
echo "1- The following command will launch the parallel workflow through $count torque jobs:"
echo qsub ${torque_script}
echo "2- The second command will merge all outputs to generate the contact maps:"
echo qsub ${torque_script_s2}
