#!/bin/bash
## Nicolas Servant
## Institut Curie

## Mapping statistics for Hi-C data

dir=$(dirname $0)

. $dir/hic.inc.sh

## Usage
function usage {
    echo -e "Usage : ./mappingstat.sh"
    echo -e "-i"" <input directory>"
    echo -e "-c"" <config>"
    echo -e "-o"" <output directory/prefix>"
    echo -e "-h"" <help>"
    exit
}

mode=global
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) ncrna_conf=$2; shift;;
	(-i) input_dir=$2; shift;;
	(-o) output_dir=$2; shift;;
	(-l) mode=local; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

read_config $ncrna_conf

make_stat()
{
    local file="$1"
    local output="$2"

    if [ ! -z "$DRY_RUN" ]; then echo "writing stats of $file to $output"; return; fi

    echo "## File $file" > ${output}
    awk '
      BEGIN {FS="\t"}
      {
        mark = $5;
        if(mark == "UN") {
          unmap++;
        } else if (mark == "R") {
          undefine++;
        } else if (mark == "U") {
          map++;
        }
      }
      
      END {
        print "Uniquely_Mapped",map;
        print "Undefined",undefine;
        print "Unmapped",unmap;
    }' $file >> ${output}
}

for r in $(get_aln_for_stats ${mode})
do
    R1=$r
    R2=$(echo $r | get_R2)
    R_STAT1=$(get_stat_file $mode $R1)
    R_STAT2=$(get_stat_file $mode $R2)
    make_stat $R1 $R_STAT1 &
    make_stat $R2 $R_STAT2 &

    wait
done
