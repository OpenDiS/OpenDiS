
.PHONY: all help clean

all: help seg_stress segseg_force loop_force 

help: 
	@echo "make sure environmental variable PYTHONPATH is present"
	@echo "you may include the following line in $$HOME/.bash_profile"
	@echo " export PYTHONPATH="
	@echo "make sure to set environmental variables correctly"
	@echo " your PYTHONPATH=$(PYTHONPATH)"
	@echo "      CTYPESGEN_DIR=$(CTYPESGEN_DIR)"

seg_stress: test_seg_stress.py
	python3 test_seg_stress.py

segseg_force: test_segseg_force.py
	python3 test_segseg_force.py

loop_force: test_loop_force.py
	python3 test_loop_force.py

clean: 
	@echo "nothing to clean"	
