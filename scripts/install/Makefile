## HiC-Pro
## Copyleft 2015,2016 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007. 

## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
VNUM = $(shell bash $(MK_PATH)/bin/HiC-Pro --version | cut -d " " -f 3)

SCRIPTS=$(MK_PATH)/scripts
SOURCES=$(SCRIPTS)/src

all : configure

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
ifndef PREFIX
ifdef prefix
PREFIX = $(prefix)
else
PREFIX = /usr/local/bin/
endif
endif

######################################
## Dependencies
##
######################################
configure: config_check
	./scripts/install/install_dependencies.sh -c $(CONFIG_SYS) -p ${PREFIX} -o $(PREFIX)/HiC-Pro_$(VNUM) -q


