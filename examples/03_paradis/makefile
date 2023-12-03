
.PHONY: all help test_paradis test_send_nodes_to_paradis clean

all: help test_paradis test_send_nodes_to_paradis

help: 
	@echo "make sure environmental variable PYTHONPATH is present"
	@echo "you may include the following line in $$HOME/.bash_profile"
	@echo " export PYTHONPATH="
	@echo "make sure to set environmental variables correctly"
	@echo " your PYTHONPATH=$(PYTHONPATH)"
	@echo "      CTYPESGEN_DIR=$(CTYPESGEN_DIR)"

test_paradis: test_paradis.py
	python3 test_paradis.py paradis_default

test_send_nodes_to_paradis: test_send_nodes_to_paradis.py
	python3 test_send_nodes_to_paradis.py paradis_default

clean: 
	rm -rf paradis_results
