#*******************************************************************************************************#
# File:              ShellScript.tmpl                                              
# Author:            hexiang                                     | Boris Jeremic,                       #
# Date:              2017-04-20 09:58:51                         | University of California, Davis,95616#
# Description:       #############                               | California                           #
# Rev:               Version 1                                   | jeremic@ucdavis.edu                  #
# Email:             hexwang@ucdavis.edu                         | Computational Geomechanics Group     #
# * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  # 
#                           Last Modified time: 2017-04-21 00:40:58                                     #              
#  * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #         
# The copyright to the computer program(s) herein is the property of Hexiang Wang and Boris Jeremic     #
# The program(s) may be used and/or copied only with written permission of Hexiang Wang or in accordance# 
# with the terms and conditions stipulated in the agreement/contract under which the program have been  #
# supplied.                                                                                             #
#*******************************************************************************************************#
TOO=g++ -std=gnu++11 -fpic -c
TOSO=g++ -shared -o 
TOOUT=g++ -std=gnu++11
BUILD_DIR = $(realpath $(CURDIR))/../build
SOURCE_DIR= $(realpath $(CURDIR))

all: rebuild

preprocess: geometry_preprocess
	@cd ..; mkdir build; 
	@chmod +x ./geometry_preprocess
	@ln -sf $(SOURCE_DIR)/geometry_preprocess $(BUILD_DIR)

node_extraction: node_extraction.py
	@chmod +x ./node_extraction.py
	@ln -sf $(SOURCE_DIR)/node_extraction.py $(BUILD_DIR)
clean:
	@rm -rf $(realpath $(CURDIR))/../build

rebuild: clean preprocess node_extraction