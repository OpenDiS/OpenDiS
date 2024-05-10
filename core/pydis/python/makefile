
# Require Ctypesgen be installed
# Install by 
#
#     git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git
#
# Add the following to ~/.bash_profile
#
#     export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git

CC_PREPROCESS   = ${CC} -E ${DEFS}

CALFORCE_HEADER_PATH = ../c/calforce
CALFORCE_HEADER_FILES = $(CALFORCE_HEADER_PATH)/SegSegForce.h $(CALFORCE_HEADER_PATH)/SegmentStress.h $(CALFORCE_HEADER_PATH)/StressDueToSeg.h $(CALFORCE_HEADER_PATH)/SegSegForce_SBN1.h $(CALFORCE_HEADER_PATH)/SegSegForce_SBN1_SBA.h

INCLUDE_HEADER_PATH = ../c/include
INCLUDE_HEADER_FILES = $(INCLUDE_HEADER_PATH)/Home.h $(INCLUDE_HEADER_PATH)/Init.h $(INCLUDE_HEADER_PATH)/ParadisProto.h

HEADER_FILES = ${CALFORCE_HEADER_FILES} ${INCLUDE_HEADER_FILES}

LIB_PYDIS_PATH  = ../../../lib
LIB_PYDIS_SO  = libpydis.so

PYDIS_LIB_PY = pydis_lib.py

CTYPESGEN_DIR = ./ctypesgen

#ifndef CTYPESGEN_DIR
#  $(info Install Ctypesgen by)
#  $(info   git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git)
#  $(info Add the following to ~/.bash_profile)
#  $(info   export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git)
#  $(error CTYPESGEN_DIR is not set)
#endif

.PHONY: all clean

all: $(PYDIS_LIB_PY)

$(PYDIS_LIB_PY): ${HEADER_FILES} ${LIB_PYDIS_PATH}/${LIB_PYDIS_SO}
	${CTYPESGEN_DIR}/run.py --cpp="${CC_PREPROCESS}" ${HEADER_FILES} -l ${LIB_PYDIS_SO} -o $@
	mv $@ ../../../lib

clean: 
	rm -rf $(PYDIS_LIB_PY) *.pyc __pycache__

