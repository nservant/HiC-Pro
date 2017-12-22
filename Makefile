## HiC-Pro
## Copyleft 2015,2016 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 Licence


## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))

INST_SCRIPTS=$(MK_PATH)/scripts
INST_SOURCES=$(INST_SCRIPTS)/src
CONFIGURE_OUT=$(wildcard ./config-system.txt)
CONFIG_SYS=$(wildcard ./config-install.txt)


install : config_check mapbuilder readstrimming iced cp

######################################
## Config file
##
######################################
config_check:
ifneq ("$(CONFIGURE_OUT)","")
include $(CONFIGURE_OUT)
else
	$(error config-system.txt file not found. Please run 'make configure' first)
endif

######################################
## Configure
##
######################################
configure:
ifneq ("$(CONFIG_SYS)","")
	make -f ./scripts/install/Makefile CONFIG_SYS=$(CONFIG_SYS)
else
	$(error config-install.txt file not found !)
endif

######################################
## Compile
##
######################################

## Build C++ code
mapbuilder: $(INST_SOURCES)/build_matrix.cpp
	(g++ -Wall -O2 -std=c++0x -o build_matrix ${INST_SOURCES}/build_matrix.cpp; mv build_matrix ${INST_SCRIPTS})

readstrimming: $(INST_SOURCES)/cutsite_trimming.cpp
	(g++ -Wall -O2 -std=c++0x -o cutsite_trimming ${INST_SOURCES}/cutsite_trimming.cpp; mv cutsite_trimming ${INST_SCRIPTS})

## Build Python lib
iced: $(INST_SOURCES)/ice_mod
	(cp $(INST_SOURCES)/ice_mod/iced/scripts/ice ${INST_SCRIPTS}; cd $(INST_SOURCES)/ice_mod/; ${PYTHON_PATH}/python setup.py install;)

test: config_check
	@echo ${PYTHON_PATH}

######################################
## Create installed version
##
######################################

cp:
ifneq ($(realpath $(MK_PATH)), $(realpath $(INSTALL_PATH)))
	cp -Ri $(MK_PATH) $(INSTALL_PATH)
endif
	@echo "HiC-Pro installed in $(realpath $(INSTALL_PATH)) !"

