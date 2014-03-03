#!/bin/bash
## Nicolas Servant
## Institut Curie

## Mapping statistics for Hi-C data

## Usage
function usage {
    echo -e "Usage : ./mappingstat.sh"
    echo -e "-i"" <input directory>"
    echo -e "-o"" <output directory/prefix>"
    echo -e "-h"" <help>"
    exit
}


################### Initialize ###################

set -- $(getopt i:o: "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-i) input_dir=$2; shift;;
	(-o) output_dir=$2; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

if [ ! -d $input_dir ]
then
    echo "$input_dir" is not a directory
    usage
    exit
fi

## U = Uniquely mapped reads
## R = Repeated
## UN = Unaligned reads

R1_STAT=${output_dir}_R1.mapstat
R2_STAT=${output_dir}_R2.mapstat

## R1 tag
echo "## Mapping statistics - input directory = $input_dir" > ${R1_STAT}
for ppfile in ${input_dir}/*R1*bowtie_final.aln
do
    echo "## File $ppfile" >> ${R1_STAT}
    awk 'BEGIN{FS="\t"}{mark=$5; if(mark=="UN"){unmap++} if(mark=="R"){undefine++} if(mark=="U"){map++}}END{print "Uniquely_Mapped",map; print "Undefined",undefine; print "Unmapped",unmap;}' $ppfile >> ${R1_STAT}
done

## R2 tag
echo "## Mapping statistics - input directory = $input_dir" > ${R2_STAT}
for ppfile in ${input_dir}/*R2*bowtie_final.aln
do
    echo "## File $ppfile" >> ${R2_STAT}
    awk 'BEGIN{FS="\t"}{mark=$5; if(mark=="UN"){unmap++} if(mark=="R"){undefine++} if(mark=="U"){map++}}END{print "Uniquely_Mapped",map; print "Undefined",undefine; print "Unmapped",unmap;}' $ppfile >> ${R2_STAT}
done