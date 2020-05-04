## hic.inc.sh
##
## Copyright (c) 2015 Institut Curie                               
## Author(s): Eric Viara, Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

###########################
## Load Configuration
###########################

set -o pipefail  # trace ERR through pipes
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

CURRENT_PATH=`dirname $0`

tmpfile1=/tmp/hic1.$$
tmpfile2=/tmp/hic2.$$
tmpmkfile=/tmp/hicmk.$$
trap "rm -f $tmpfile1 $tmpfile2 $tmpmkfile" 0 1 2 3

abspath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

filter_config()
{
    sed -e 's/#.*//' | egrep '^[ \t]*[a-zA-Z_][a-zA-Z0-9_]*[ \t]*:?=' | sed -e 's/[ \t]*:=[ \t]*/ :=/' -e 's/[ \t][^:]*=[ \t]*/ =/' -e 's/\([^ \t]*\)=/\1 =/' -e 's/ *$//g' | sort -u -k 1b,1
}

 read_config()
{
    local conf=$1
    cat $conf > $tmpmkfile
    echo "_dummy_target_:" >> $tmpmkfile
    make -f $tmpmkfile -p -n | filter_config > $tmpfile1
    cat $conf | filter_config > $tmpfile2

    eval "$(join $tmpfile1 $tmpfile2 | awk -F' =' '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)}')"

    ## Define BOWTIE outputs
    BOWTIE2_IDX=${BOWTIE2_IDX_PATH}/${REFERENCE_GENOME}; export BOWTIE2_IDX
    BOWTIE2_GLOBAL_OUTPUT_DIR=$BOWTIE2_OUTPUT_DIR/bwt2_global; export BOWTIE2_GLOBAL_OUTPUT_DIR
    BOWTIE2_LOCAL_OUTPUT_DIR=$BOWTIE2_OUTPUT_DIR/bwt2_local; export BOWTIE2_LOCAL_OUTPUT_DIR
    BOWTIE2_FINAL_OUTPUT_DIR=$BOWTIE2_OUTPUT_DIR/bwt2; export BOWTIE2_FINAL_OUTPUT_DIR
    
    ## Clean RAW_DIR variable
    RAW_DIR=$(echo $RAW_DIR | sed -e 's|^\./||')
}

## Load System config
SYS_CONF=$CURRENT_PATH/../config-system.txt
if [ -e "$SYS_CONF" ]; then
    read_config $SYS_CONF
else
    echo "Error - System config file not found"
    exit
fi

## load Hi-C config
if [ ! -z "$CONF" ]; then
    CONF=`abspath $CONF`
    if [ -e "$CONF" ]; then
	read_config $CONF
    else
	echo "Error - Hi-C config file '$CONF' not found"
	exit
    fi
fi

###########################
## Subroutine for scripts
###########################

die() 
{ 
    echo "Exit: $@" 1>&2 
    exit 1
}

exec_cmd()
{
    echo $*
    if [ -z "$DRY_RUN" ]; then
	eval "$@" 
    fi
}

exec_ret()
{
    if [ -z "$DRY_RUN" ]; then
	eval "$@" 
    fi
}

add_ext()
{
    local file=$1
    local ext=$2
    echo $file | grep "\${ext}$" > /dev/null
    if [ $? = 0 ]; then
	echo ${file}
    else
	echo ${file}${ext}
    fi
}

add_fastq()
{
    add_ext $1 $2
}

get_R1()
{
    sed -e "s/${PAIR2_EXT}/${PAIR1_EXT}/"
}

get_R2()
{
    sed -e "s/${PAIR1_EXT}/${PAIR2_EXT}/"
}

get_pairs()
{
    sed -e "s/${PAIR1_EXT}//;s/${PAIR2_EXT}//"
}

filter_pairs()
{
	get_R1 | sort -u
}

get_data_type()
{
    ## return the highest possible input files type
    nb_fq=$(find -L $RAW_DIR -mindepth 2 -maxdepth 2 -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz"| wc -l)
    nb_bam=$(find -L $RAW_DIR -mindepth 2 -maxdepth 2 -name "*.bam" -o -name "*.sam" | wc -l)
    nb_vpairs=$(find -L $RAW_DIR -mindepth 2 -maxdepth 2 -name "*.validPairs" | wc -l)
    nb_allvpairs=$(find -L $RAW_DIR -mindepth 2 -maxdepth 2 -name "*.allValidPairs" | wc -l)
    nb_mat=$(find -L $RAW_DIR -mindepth 2 -maxdepth 4 -name "*.matrix" | wc -l)

    if (( $nb_mat > 0 )); then
        INPUT_DATA_TYPE="mat"
    elif (( $nb_allvpairs > 0 )); then
        INPUT_DATA_TYPE="allvalid"
    elif (( $nb_vpairs > 0 )); then
        INPUT_DATA_TYPE="valid"
    elif (( $nb_bam > 0 )); then
        INPUT_DATA_TYPE="bam"
    elif (( $nb_fq > 0 )); then
        INPUT_DATA_TYPE="fastq"
    else
	die "Error in input type.'.fastq|.fq|.bam|.validPairs|.allValidPairs|.matrix' files are expected." #!
    fi
    echo $INPUT_DATA_TYPE
}

#
# function called by get_hic_files, using "local" variables of get_hic_files
#

set_ext2fastq()
{
    local file=$1
    local ext=$2
    file=$(echo $file | sed -e "s/\.fastq$//" -e "s/\.fastq.gz$//" -e "s/\.fq$//" -e "s/\.fq.gz$//") #!
    echo ${file}${ext}
}

get_hic_files_build_list()
{
    local file=$idir/$(set_ext2fastq $fastq $ext)
    if [ ! -z "$list" ]; then
	list="$list
$file"
    else
	list=$file
    fi
}

filter_rawdir()
{
    sed -e "s|\./${RAW_DIR}/||g" -e "s|${RAW_DIR}/||g"
}

get_hic_files()
{
    local idir=$1
    local ext=$2
    if [ ! -z "$PBS_ARRAYID" ]; then TASKID=$PBS_ARRAYID; fi
    if [ ! -z "$SGE_TASK_ID" ]; then TASKID=$SGE_TASK_ID; fi
    if [ ! -z "$SLURM_ARRAY_TASK_ID" ]; then TASKID=$SLURM_ARRAY_TASK_ID; fi
    if [ ! -z "$LSB_JOBINDEX" ]; then TASKID=$LSB_JOBINDEX; fi
    if [ ! -z "$FASTQFILE" ]; then
	if [ ! -z "$TASKID" ]; then
	    local input_data_type=$(get_data_type)
	    ## deal with fq/fastq extension
	    if [ ${input_data_type} == "fastq" ]; then
		pattern=".fastq(.gz)*$|.fq(.gz)*$"
	    else
		pattern=".${input_data_type}$"
	    fi
	    ## raw data for mapping
	    if [[ $ext == ".fastq" || $ext == ".fq" ]]; then
                cat $FASTQFILE | filter_rawdir | filter_pairs | awk "NR == $TASKID && \$1 ~ \"${ext}(.gz)*$\"{printf(\"%s/%s${ext}\n\", \"$idir\", gensub(\"${ext}(.gz)*$\", \"\", \$1));}"
	    else
		cat $FASTQFILE | filter_rawdir | filter_pairs | awk "NR == $TASKID {printf(\"%s/%s${ext}\n\", \"$idir\", gensub(\"${pattern}\", \"\", \$1));}"
    	    fi
	    return
	fi
	local list=
	for fastq in $(cat $FASTQFILE | filter_rawdir ); do
	    if [[ ${ext} == ".fastq" || ${ext} == ".fq" ]]
	    then
		if [[ $fastq =~ "${ext}" ]]
		then
		    get_hic_files_build_list
		fi
	    else
		get_hic_files_build_list
	    fi
	done
	echo "$list" | filter_pairs
    elif [ ! -z "$FASTQLIST" ]; then
	local list=
	for fastq in $(echo $FASTQLIST | filter_rawdir | sed -e 's/[,;]/ /g'); do
           if [[ ${ext} == ".fastq" || ${ext} == ".fq" ]]
           then
               if [[ $fastq =~ "${ext}" ]]
               then
                   get_hic_files_build_list
               fi
           else
               get_hic_files_build_list
           fi
	done
	echo "$list" | filter_pairs
    else
	find $idir \( -name \*${ext} -o -name \*${ext}.gz \) -follow -print | filter_pairs
    fi
}

get_sample_dir()
{
    local file=$1
    echo $(basename $(dirname ${file}))
}

get_fastq_for_bowtie_global()
{
    local input_data_type=$(get_data_type)    
    if [[ $input_data_type == "fastq" ]]
    then
        ifastq=$(get_hic_files $RAW_DIR .fastq | grep "$PAIR1_EXT")
        ifq=$(get_hic_files $RAW_DIR .fq | grep "$PAIR1_EXT")
	echo "$ifastq $ifq"
    fi
}

get_fastq_for_bowtie_local()
{
    get_hic_files $BOWTIE2_GLOBAL_OUTPUT_DIR _${REFERENCE_GENOME}.bwt2glob.unmap.fastq
}

get_bam_for_pp()
{
    local=$1
    if [ "$local" = 1 ]; then
	get_hic_files $BOWTIE2_LOCAL_OUTPUT_DIR _${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam
    else
	get_hic_files $BOWTIE2_GLOBAL_OUTPUT_DIR _${REFERENCE_GENOME}.bwt2glob.bam
    fi
}

get_global_aln_for_stats()
{
    get_hic_files ${BOWTIE2_GLOBAL_OUTPUT_DIR} _${REFERENCE_GENOME}.bwt2glob.bam
}

get_local_aln_for_stats()
{
    get_hic_files ${BOWTIE2_LOCAL_OUTPUT_DIR} _${REFERENCE_GENOME}.bwt2glob.unmap_bwt2loc.bam
}

get_stat_file()
{
    local file=$1
    local sample_dir=$(get_sample_dir ${file})
    local prefix=$(echo ${sample_dir}/$(basename $file) | sed -e 's/.bwt2glob.bam//')

    echo ${BOWTIE2_FINAL_OUTPUT_DIR}/$prefix.mapstat
}

get_bam_from_raw_dir()
{
    get_hic_files $RAW_DIR .bam | grep "$PAIR1_EXT"
}

get_sam_for_merge()
{
    local input_data_type=$(get_data_type)
    if [[ $input_data_type == "fastq" ]]
    then
	bam=$(get_hic_files ${BOWTIE2_FINAL_OUTPUT_DIR} _${REFERENCE_GENOME}.bwt2merged.bam)
    elif [[ $input_data_type == "bam" ]]
    then
	bam=$(get_bam_from_raw_dir)
    fi
    echo $bam
}

get_sam_for_combine()
{
    get_hic_files ${BOWTIE2_GLOBAL_OUTPUT_DIR} _${REFERENCE_GENOME}.bwt2glob.bam   
}

get_paired_bam()
{
    get_hic_files ${BOWTIE2_FINAL_OUTPUT_DIR} _${REFERENCE_GENOME}.bwt2pairs.bam | get_R1 | sed -e "s/${PAIR1_EXT}//" -e "s/_${REFERENCE_GENOME}.bwt2merged//"
}

