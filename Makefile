## HiC-Pro
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007. 

## DO NOT EDIT THE REST OF THIS FILE!!

SCRIPTS=./scripts
SOURCES=$(SCRIPTS)/src

all : install

install : checkdep mapbuilder readstrimming iced

######################################
## Config file
##
######################################
config_check:
ifndef CONFIG_SYS
	$(error CONFIG_SYS is not defined. Please run 'make CONFIG_SYS=config-install.txt install')
else		
include $(CONFIG_SYS)
endif

######################################
## Dependencies
##
######################################
checkdep: config_check
	./scripts/install/install_dependencies.sh -c $(CONFIG_SYS)

######################################
## Compile
##
######################################

## Build C++ code
mapbuilder: $(SOURCES)/build_matrix.cpp
	(g++ -Wall -O2 -std=c++0x -o build_matrix ${SOURCES}/build_matrix.cpp; mv build_matrix ${SCRIPTS})

readstrimming: $(SOURCES)/cutsite_trimming.cpp
	(g++ -Wall -O2 -std=c++0x -o cutsite_trimming ${SOURCES}/cutsite_trimming.cpp; mv cutsite_trimming ${SCRIPTS})

## Build Python lib
iced: $(SOURCES)/ice_mod
	(cp $(SOURCES)/ice_mod/iced/scripts/ice ${SCRIPTS}; cd $(SOURCES)/ice_mod/; python setup.py install --user;)


