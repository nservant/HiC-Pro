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
    echo -e "-o"" <installation folder>"
    echo -e "-q"" <quiet>"
    echo -e "-h"" <help>"
    exit;
}

echo -e "$RED""Make sure internet connection works for your shell prompt under current user's privilege ...""$NORMAL";
echo -e "$BLUE""Starting $SOFT installation ...""$NORMAL";


################### Initialize ###################
quiet=0
set -- $(getopt c:o:qh "$@")
while [ $# -gt 0 ]
do
    case "$1" in
	(-c) conf=$2; shift;;
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
    vercomp $pver "2.7.0"
    if [[ $? == 2 ]]; then
	echo -e "$RED""Python v2.7.0 or higher is needed [$pver detected].""$NORMAL"
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

PREFIX_BIN=/usr/local/bin

if [ ! -w $PREFIX_BIN ]; then
    PREFIX_BIN=${HOME}/bin;
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
echo "Checking dependencies ... "

wasInstalled=0;
echo "Checking Python libraries ..."
python ../scripts/install/check_pythonlib.py > install_packages_check.Rout
if [ $? == "0" ]; then
    echo -e "$BLUE""The required Python libraries appear to be already installed. ""$NORMAL"
    wasInstalled=1;
else
    echo -e "$RED""Can not proceed without the required Python libraries, please install them and re-run""$NORMAL"
    exit 1;
fi


################  R  ###################

wasInstalled=0;
which R > /dev/null;
if [ $? == "0" ]; then
    echo "Checking R installation ..."
    R CMD BATCH ../scripts/install/check_Rpackages.R > check_Rpackages.Rout
    check=`grep proc.time check_Rpackages.Rout`;
    if [ $? == "0" ]; then
	echo -e "$BLUE""The required R packages appear to be already installed. ""$NORMAL"
	wasInstalled=1;
    fi
else
    echo -e "$RED""Can not proceed without R, please install and re-run""$NORMAL"
    exit 1;
fi

#Install R Packages
if [ $wasInstalled == 0 ]; then
    echo "Installing missing R packages ..."
    R CMD BATCH ../scripts/install/install_Rpackages.R install_Rpackages.Rout

    R CMD BATCH ../scripts/install/check_Rpackages.R > check_Rpackages.Rout
    check=`grep proc.time check_Rpackages.Rout`;
    if [ $? == "0" ]; then
	echo -e "$BLUE""R packages appear to be installed successfully""$NORMAL"
    else
	echo -e "$RED""R packages NOT installed successfully. Look at the tmp/install_Rpackages.Rout for additional informations""$NORMAL"; exit 1;
    fi
fi

################ Bowtie2 ###################

wasInstalled=0;
which bowtie2 > /dev/null;
if [ $? = "0" ]; then
	echo -e "$BLUE""Bowtie2 Aligner appears to be already installed. ""$NORMAL"
	wasInstalled=1;
fi


if [ $wasInstalled == 0 ]; then
    echo "Installing Bowtie2 ..."
    $get bowtie2-2.2.4-source.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-source.zip/download
    unzip bowtie2-2.2.4-source.zip
    cd bowtie2-2.2.4
    make
    cd ..
    mv bowtie2-2.2.4 $PREFIX_BIN
    export PATH=$PREFIX_BIN/bowtie2-2.2.4/:$PATH
    wasInstalled=0;
fi
 
#some bowtie tests
if [ $wasInstalled == 0 ]; then
    check=`bowtie2 --version 2>&1`;
    if [ $? = "0" ]; then
	echo -e "$BLUE""Bowtie2 Aligner appears to be installed successfully""$NORMAL"
    else
	echo -e "$RED""Bowtie2 Aligner NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

################ samtools  ###################

wasInstalled=0;
which samtools > /dev/null
if [ $? = "0" ]; then
        
    samver=`samtools 2>&1 | grep Version | cut -d" " -f2`
    vercomp $samver "1.0"
    if [[ $? == 2 ]]; then
	echo -e "$RED""samtools v1.0 or higher is needed [$samver detected].""NORMAL"
	exit 1;
    fi

	echo -e "$BLUE""Samtools appears to be already installed. ""$NORMAL"
	wasInstalled=1;
fi

if [ $wasInstalled == 0 ]; then
    echo "Installing samtools ..."
    #From sources
    $get samtools-1.1.tar.bz2  http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2/download
    tar -xvjpf samtools-1.1.tar.bz2
    cd samtools-1.1
    make
    cd ..
    mv samtools-1.1 $PREFIX_BIN
    export PATH=$PREFIX_BIN/samtools-1.1/:$PATH
    wasInstalled=0;
fi

if [ $wasInstalled == 0 ]; then
    check=`samtools view -h 2>&1 | grep -i options`;
    if [ $? = "0" ]; then
	echo -e "$BLUE""samtools appears to be installed successfully""$NORMAL"
    else
	echo -e "$RED""samtools NOT installed successfully""$NORMAL"; exit 1;
    fi
fi

## Clean up
cd ..
rm -rf ./tmp

echo -e "$RED""Dependencies checked !""$NORMAL"

################ Create the config-system file ###################
CUR_DIR=`pwd`
echo -e "$BLUE""Check $SOFT configuration ... ""$NORMAL"

echo "#######################################################################" > config-system.txt
echo "## $SOFT - System settings" >> config-system.txt
echo "#######################################################################" >> config-system.txt

echo "#######################################################################" >> config-system.txt
echo "## Required Software - Specified the DIRECTORY name of the executables" >> config-system.txt
echo "## If not specified, the program will try to locate the executable" >> config-system.txt
echo "## using the 'which' command" >> config-system.txt
echo "#######################################################################" >> config-system.txt

which R > /dev/null
if [ $? = "0" ]; then
    echo "R_PATH = "`dirname $(which R)` >> config-system.txt
else
    die "R_PATH not found. Exit." 
fi

which bowtie2 > /dev/null
if [ $? = "0" ]; then
    echo "BOWTIE2_PATH = "`dirname $(which bowtie2)`  >> config-system.txt
else
    die "BOWTIE2_PATH not found. Exit." 
fi

which samtools > /dev/null
if [ $? = "0" ]; then
    echo "SAMTOOLS_PATH = "`dirname $(which samtools)`  >> config-system.txt
else
    die "SAMTOOLS_PATH not found. Exit." 
fi

which python > /dev/null
if [ $? = "0" ]; then
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
    echo -e "$RED""Warning : Scheduler system not defined - Default is Torque/PBS""$NORMAL"
    CLUSTER_SYS="TORQUE"; 
fi
if [ $CLUSTER_SYS == "TORQUE" ]; then 
    #ln -s scripts/make_torque_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_torque_script.sh" >> config-system.txt
    echo -e "$BLUE""Configuration for TORQUE/PBS system.""$NORMAL"
elif [ $CLUSTER_SYS == "SGE" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_sge_script.sh" >> config-system.txt
    echo -e "$BLUE""Configuration for SGE system.""$NORMAL"
elif [ $CLUSTER_SYS == "SLURM" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_slurm_script.sh" >> config-system.txt
    echo -e "$BLUE""Configuration for SLURM system.""$NORMAL"
elif [ $CLUSTER_SYS == "LSF" ]; then 
    #ln -s scripts/make_sge_scripts.sh scripts/make_cluster_scripts.sh
    echo "CLUSTER_SCRIPT = ${install_dir}/scripts/make_lsf_script.sh" >> config-system.txt
    echo -e "$BLUE""Configuration for LSF system.""$NORMAL"
else
    die "$CLUSTER_SYS unknown. 'TORQUE/SGE/SLURM/LSF' systems are currently supported. Please change the CLUSTER_SYS variable and re-run the installation process. Exit."
fi

## check rights in PREFIX folder
if [[ -z $PREFIX ]]; then PREFIX=/usr/local/bin; fi
if [ ! -w $PREFIX ]; then
    die "Cannot install HiCPro in $PREFIX directory. Maybe missing super-user (root) permissions to write there. Please specify another directory in the config-install.txt file (PREFIX=)";
fi 

echo ;
## End of dependencies check ##
