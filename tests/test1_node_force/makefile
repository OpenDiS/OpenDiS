
.PHONY: all help run clean

all: help run

# modify the PYTHONPATH variable
PYTHONPATH=$(abspath $(CURDIR)/../../lib):$(abspath $(CURDIR)/../../python):

help: 
	@echo "make sure environmental variable PYTHONPATH is present"
	@echo "you may include the following line in $$HOME/.bash_profile"
	@echo " export PYTHONPATH="
	@echo "make sure to set environmental variables correctly"
	@echo " your PYTHONPATH=$(PYTHONPATH)"
	@echo "      CTYPESGEN_DIR=$(CTYPESGEN_DIR)"

segseg_force: test_segseg_force.py
	python3 test_segseg_force.py

clean: 
	@echo "nothing to clean"	
