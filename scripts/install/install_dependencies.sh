#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

## This script aims in installing most of the dependies of the HiC-Pro tool.
## Serval checks are done to ensure compilation of code.
##


NORMAL="\\033[0;39m"
RED="\\033[0;31m"
BLUE="\\033[0;34m"
GREEN="\\033[0;32m"
YELLOW="\\033[1;33m"
SOFT="HiC-Pro"

## 0 =
## 1 >
## 2 <
vercomp () {
    if [[ $1 == $2 ]]
    then
        return 0
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done

    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [[ -z ${ver2[i]} ]]
        then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if ((10#${ver1[i]} > 10#${ver2[i]}))
        then
            return 1
        fi
        if ((10#${ver1[i]} < 10#${ver2[i]}))
        then
            return 2
        fi
    done
    return 0
}


die() {
    echo -e "$RED""Exit - ""$*""$NORMAL" 1>&2
    exit 1
}

function usage {
    echo -e "Usage : ./install_all.sh"
    echo -e "-c"" <configuration install file>"
    echo -e "-p"" <prefix>"
    echo -e "-o"" <installation folder>"
    echo -e "-q"" <quiet>"
    echo -e "-h"" <help>"
    exit;
}

echo -e "$YELLOW""Make sure internet connection works for your shell prompt under current user's privilege ...""$NORMAL";
echo -e "$NORMAL""Starting $SOFT installation !""$NORMAL";


################### Initialize ###################
quiet=0
set -- $(getopt c:p:o:qh "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf=$2; shift;;
	(-p) prefix=$2 shift;;
	(-o) install_dir=$2; shift;;
	(-q) quiet=1; shift;;
	(-h) usage;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
	(*)  break;;
    esac
    shift
done

if [[ -z $install_dir ]] 
then 
    die "Error : Installation directory not defined (-o)"
fi
if [[ ! -e $conf ]]
then 
    die "Error : Configuration file not found"
fi


################### Read the config file ###################

while read curline_read; do
    curline=${curline_read// /}
    if [[ $curline != \#* && ! -z $curline ]]; then
	var=`echo $curline | awk -F= '{print $1}'`
	val=`echo $curline | awk -F= '{print $2}'`

	if [[ $var =~ "_PATH" ]]
	then
	    if [[ ! -z $val ]]; then
		echo "export $val in PATH"
		export PATH=$val:$PATH
	    fi
	else
	    export $var=$val
	fi
    fi
done < $conf

################### Search standard tools ###################

#check for make
which make > /dev/null;
if [ $? != "0" ]; then
	echo -e "$RED""Can not proceed without make, please install and re-run (Mac users see: http://developer.apple.com/technologies/xcode.html)""$NORMAL"
	exit 1;
fi

#check for g++
which g++ > /dev/null;
if [ $? != "0" ]; then
	echo -e "$RED""Can not proceed without g++, please install and re-run""$NORMAL"
	exit 1;
fi

# check for unzip (bowtie)
which unzip > /dev/null;
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without unzip, please install and re-run""$NORMAL"
    exit 1;
fi

# python
which python > /dev/null;
if [ $? != "0" ]; then
    echo -e "$RED""Can not proceed without Python, please install and re-run""$NORMAL"
    exit 1;
else
    pver=`python --version 2>&1 | cut -d" " -f2`
    vercomp $pver "3.7.0"
    if [[ $? == 2 ]]; then
	echo -e "$RED""Python v3.7.X or higher is needed [$pver detected].""$NORMAL"
	exit 1;
    fi
fi

 

#check OS (Unix/Linux or Mac)
os=`uname`;

# get the right download program
if [ "$os" = "Darwin" ]; then
	# use curl as the download program 
	get="curl -L -o"
else
	# use wget as the download program
	get="wget --no-check-certificate -O"
fi

if [ -d ./tmp ]; then
    rm -r ./tmp
fi
mkdir ./tmp
cd ./tmp

################ Install dependencies  ###################

## By default, dependencies will be installed in the same path than HiC-Pro
if [[ ! -z ${prefix} ]]; then
    PREFIX=${prefix}
fi
PREFIX_BIN=${PREFIX}

if [ ! -w $PREFIX_BIN ]; then
   quiet=1
fi
if [[ $quiet == 0 ]]; then
    echo "Where should missing software prerequisites be installed ? [$PREFIX_BIN] "
    read ans
    ans=${ans:-$PREFIX_BIN}
    PREFIX_BIN=$ans
fi

if [ ! -d $PREFIX_BIN ]; then
    echo "Directory $PREFIX_BIN does not exist!"

    if [[ $quiet == 0 ]]; then
	echo -n "Do you want to create $PREFIX_BIN folder ? (y/n) [n] : "
	read ans
	if [ XX${ans} = XXy ]; then
            mkdir $PREFIX_BIN || die "Cannot create  $PREFIX_BIN folder. Maybe missing super-user (root) permissions"
	else
            die "Must specify a directory to install required softwares!"
	fi
    else
	die "Error - unable to install/check dependancies !"
    fi
fi

if [ ! -w $PREFIX_BIN ]; then
    die "Cannot write to directory $PREFIX_BIN. Maybe missing super-user (root) permissions to write there.";
fi 

################  Python  ###################
echo "Checking dependencies"

wasInstalled=0;
echo -n "- Python libraries ..."
python ../scripts/install/check_pythonlib.py > check_python.log 2>&1
if [ $? == "0" ]; then
    echo -e "$GREEN""OK""$NORMAL"
    wasInstalled=1;
else
    echo -e "$RED""\nCan not proceed without the required Python libraries, please install them and re-run""$NORMAL"
    exit 1;
fi


################  R  ###################
wasInstalled=0;
echo -n "- R installation ..."
which R > /dev/null;
if [ $? == "0" ]; then
    R CMD BATCH ../scripts/install/check_Rpackages.R > check_Rpackages.log
    check=`grep proc.time check_Rpackages.Rout`;
    if [ $? == "0" ]; then
	echo -e "$GREEN""OK""$NORMAL"
	wasInstalled=1;
    fi
else
    echo -e "$RED""\nCan not proceed without R, please install and re-run""$NORMAL"
    exit 1;
fi

#Install R Packages
if [ $wasInstalled == 0 ]; then
    echo -n "\n  -- Installing missing R packages ..."
    R CMD BATCH ../scripts/install/install_Rpackages.R install_Rpackages.Rout
    R CMD BATCH ../scripts/install/check_Rpackages.R > check_Rpackages.Rout
    check=`grep proc.time check_Rpackages.Rout`;
    if [ $? == "0" ]; then
	echo -e "$GREEN""OK""$NORMAL"
    else
	echo -e "$RED""\nR packages NOT installed successfully. Look at the tmp/install_Rpackages.Rout for additional informations""$NORMAL"; exit 1;
    fi
fi

################ Bowtie2 ###################
wasInstalled=0;
echo -n "- Bowtie2 installation ..."
which bowtie2 > /dev/null;
if [ $? == "0" ]; then
    echo -e "$GREEN""OK""$NORMAL"
    wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo -e -n "\n -- Installing Bowtie2 ..."
    $get bowtie2-2.3.5-source.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5/bowtie2-2.3.5-source.zip/download
    unzip bowtie2-2.3.5-source.zip
    cd bowtie2-2.3.5
    make
    cd ..
    mv bowtie2-2.3.5 $PREFIX_BIN
    export PATH=$PREFIX_BIN/bowtie2-2.3.5/:$PATH
    wasInstalled=0;
fi
 
#some bowtie tests
if [ $wasInstalled == 0 ]; then
    check=`bowtie2 --version 2>&1`;
    if [ $? = "0" ]; then
	echo -e "$GREEN""OK""$NORMAL"
    else
	echo -e "$RED""\nBowtie2 Aligner NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

################ samtools  ###################

wasInstalled=0;
echo -n "- Samtools installation ..."
which samtools > /dev/null
if [ $? == "0" ]; then
    samver=`samtools 2>&1 | grep Version | cut -d" " -f2`
    vercomp $samver "1.0"
    if [[ $? == 2 ]]; then
	echo -e "$RED""\nsamtools v1.0 or higher is required [$samver detected].""NORMAL"
	exit 1;
    fi
    echo -e "$GREEN""OK""$NORMAL"
    wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo -n "\n  -- Installing samtools ..."
    #From sources
    $get samtools-1.10.tar.bz2  http://sourceforge.net/projects/samtools/files/samtools/1.10/samtools-1.10.tar.bz2/download
    tar -xvjpf samtools-1.10.tar.bz2
    cd samtools-1.10
    make
    cd ..
    mv samtools-1.10 $PREFIX_BIN
    export PATH=$PREFIX_BIN/samtools-1.10/:$PATH
    wasInstalled=0;
fi

if [ $wasInstalled == 0 ]; then
    check=`samtools view -h 2>&1 | grep -i options`;
    if [ $? = "0" ]; then
	echo -e "$GREEN""OK""$NORMAL"
    else
	echo -e "$RED""\nsamtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

## Clean up
cd ..
rm -rf ./tmp

#echo -e "$GREEN""\nDependencies checked !""$NORMAL"

################ Create the config-system file ###################
CUR_DIR=`pwd`
echo -e "$NORMAL""\nChecking $SOFT configuration""$NORMAL"

echo "#######################################################################" > config-system.txt
echo "## $SOFT - System settings" >> config-system.txt
echo "#######################################################################" >> config-system.txt

echo "#######################################################################" >> config-system.txt
echo "## Required Software - Specified the DIRECTORY name of the executables" >> config-system.txt
echo "## If not specified, the program will try to locate the executable" >> config-system.txt
echo "## using the 'which' command" >> config-system.txt
echo "#######################################################################" >> config-system.txt

which R > /dev/null
if [ $? == "0" ]; then
    echo "R_PATH = "`dirname $(which R)` >> config-system.txt
else
    die "R_PATH not found. Exit." 
fi

which bowtie2 > /dev/null
if [ $? == "0" ]; then
    echo "BOWTIE2_PATH = "`dirname $(which bowtie2)`  >> config-system.txt
else
    die "BOWTIE2_PATH not found. Exit." 
fi

which samtools > /dev/null
if [ $? == "0" ]; then
    echo "SAMTOOLS_PATH = "`dirname $(which samtools)`  >> config-system.txt
else
    die "SAMTOOLS_PATH not found. Exit." 
fi

which python > /dev/null
if [ $? == "0" ]; then
    echo "PYTHON_PATH = "`dirname $(which python)` >> config-system.txt
else
    die "PYTHON_PATH not found. Exit."
fi

echo "INSTALL_PATH = ${install_dir}" >> config-system.txt
echo "SCRIPTS = ${install_dir}/scripts" >> config-system.txt
echo "SOURCES = ${install_dir}/scripts/src" >> config-system.txt
echo "ANNOT_DIR = ${install_dir}/annotation" >> config-system.txt

## deal with scheduler system
if [ -z "$CLUSTER_SYS" ]; then 
    echo -e "$YELLOW""Warning : Scheduler system not defined - Default is Torque/PBS""$NORMAL"
    CLUSTER_SYS="TORQUE"; 
fi
if [ $CLUSTER_SYS == "TORQUE" ]; then 
    #ln -s scripts/make_torque_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_torque_script.sh" >> config-system.txt
    echo -n "- Configuration for TORQUE/PBS system ..."
elif [ $CLUSTER_SYS == "SGE" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_sge_script.sh" >> config-system.txt
    echo -n "- Configuration for SGE system ..."
elif [ $CLUSTER_SYS == "SLURM" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_slurm_script.sh" >> config-system.txt
    echo -n "- Configuration for SLURM system ..."
elif [ $CLUSTER_SYS == "LSF" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_lsf_script.sh" >> config-system.txt
    echo -n "- Configuration for LSF system ..."
else
    die "$CLUSTER_SYS unknown. 'TORQUE/SGE/SLURM/LSF' systems are currently supported. Please change the CLUSTER_SYS variable and re-run the installation process. Exit."
fi

## check rights in PREFIX folder
#if [[ -z $PREFIX ]]; then PREFIX=/usr/local/bin; fi
if [ ! -w $PREFIX ]; then
    die "Cannot install HiCPro in $PREFIX directory. Maybe missing super-user (root) permissions to write there. Please specify another directory in the config-install.txt file (PREFIX=)";
fi 
echo -e "$GREEN""OK""$NORMAL"
echo -e "$GREEN""\ndone !""$NORMAL"
## End of dependencies check ##
