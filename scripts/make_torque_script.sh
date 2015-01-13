#!/bin/bash
## Eric Viara updated 2014-05-05
##

dir=$(dirname $0)

##. $dir/hic.inc.sh

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

torque_script=HiC_torque_${PBS_SUFFIX}.sh
PPN=$(( ${N_CPU} * 2))
cat > ${torque_script} <<EOF
#!/bin/bash
#PBS -l nodes=1:ppn=${PPN},mem=${PBS_MEM},walltime=${PBS_WALLTIME}
#PBS -M ${PBS_MAIL}
#PBS -m ae
#PBS -j eo
#PBS -N HiCpro_${PBS_SUFFIX}
#PBS -q ${PBS_QUEUE}
#PBS -V
#PBS -t 1-$count

cd \$PBS_O_WORKDIR

FASTQFILE=\$PBS_O_WORKDIR/$fastqfile; export FASTQFILE
make CONFIG_FILE=${conf_file} all_qsub
EOF

chmod +x ${torque_script}

echo "The following command will launch $count torque jobs:"
echo qsub ${torque_script}
