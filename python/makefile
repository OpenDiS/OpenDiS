
# Require Ctypesgen be installed
# Install by 
#
#     git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git
#
# Add the following to ~/.bash_profile
#
#     export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git

CC_PREPROCESS   = ${CC} -E ${DEFS}

HEADER_FILES = ../src/SegSegForce.h ../src/SegmentStress.h ../src/StressDueToSeg.h ../src/SegSegForce_SBN1.h ../src/SegSegForce_SBN1_SBA.h
LIB_OPENDIS  = ../lib/libopendis.so

# Add _gpu to target name
OPENDIS_PY = opendis.py

ifndef CTYPESGEN_DIR
  $(info Install Ctypesgen by)
  $(info   git clone https://github.com/davidjamesca/ctypesgen.git ${HOME}/Codes/Ctypesgen.git)
  $(info Add the following to ~/.bash_profile)
  $(info   export CTYPESGEN_DIR=${HOME}/Codes/Ctypesgen.git)
  $(error CTYPESGEN_DIR is not set)
endif

.PHONY: all clean

all: $(OPENDIS_PY)

$(OPENDIS_PY): ${HEADER_FILES} ${LIB_OPENDIS}
	${CTYPESGEN_DIR}/run.py --cpp="${CC_PREPROCESS}" ${HEADER_FILES} -l ${LIB_OPENDIS} -o $@
	mv $@ ../lib

clean: 
	rm -rf $(OPENDIS_PY) *.pyc __pycache__

