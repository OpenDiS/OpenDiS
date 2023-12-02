
# Require Ctypesgen be installed
# Install by 
#
#     git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git
#
# Add the following to ~/.bash_profile
#
#     export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git

CC_PREPROCESS   = ${CC} -E ${DEFS}

HEADER_PATH = ../core/pydis/calforce
HEADER_FILES = $(HEADER_PATH)/SegSegForce.h $(HEADER_PATH)/SegmentStress.h $(HEADER_PATH)/StressDueToSeg.h $(HEADER_PATH)/SegSegForce_SBN1.h $(HEADER_PATH)/SegSegForce_SBN1_SBA.h
LIB_PYDIS_SO  = ../lib/libpydis.so

PYDIS_LIB_PY = pydis_lib.py

ifndef CTYPESGEN_DIR
  $(info Install Ctypesgen by)
  $(info   git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git)
  $(info Add the following to ~/.bash_profile)
  $(info   export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git)
  $(error CTYPESGEN_DIR is not set)
endif

.PHONY: all clean

all: $(PYDIS_LIB_PY)

$(PYDIS_LIB_PY): ${HEADER_FILES} ${LIB_PYDIS_SO}
	${CTYPESGEN_DIR}/run.py --cpp="${CC_PREPROCESS}" ${HEADER_FILES} -l ${LIB_PYDIS_SO} -o $@
	mv $@ ../lib

clean: 
	rm -rf $(PYDIS_LIB_PY) *.pyc __pycache__

