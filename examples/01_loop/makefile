
.PHONY: all help run clean

all: help run

help: 
	@echo "make sure environmental variable PYTHONPATH is present"
	@echo "you may include the following line in $$HOME/.bash_profile"
	@echo " export PYTHONPATH="
	@echo "make sure to set environmental variables correctly"
	@echo " your PYTHONPATH=$(PYTHONPATH)"
	@echo "      CTYPESGEN_DIR=$(CTYPESGEN_DIR)"

run: test_disnet_loop.py
	python3 -i test_disnet_loop.py

clean: 
	rm -rf test_loop_results
