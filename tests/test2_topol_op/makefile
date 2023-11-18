
.PHONY: all help clean

all: help loop_topol_op

help: 
	@echo "make sure environmental variable PYTHONPATH is present"
	@echo "you may include the following line in $$HOME/.bash_profile"
	@echo " export PYTHONPATH="
	@echo "make sure to set environmental variables correctly"
	@echo " your PYTHONPATH=$(PYTHONPATH)"
	@echo "      CTYPESGEN_DIR=$(CTYPESGEN_DIR)"

loop_topol_op: test_loop_topol_op.py
	python3 test_loop_topol_op.py

clean: 
	@echo "nothing to clean"	
