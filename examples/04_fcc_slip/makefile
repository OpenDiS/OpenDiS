# Run the tests in this folder.
#
#   make                         every test except the slow leak measurement
#   make leak                    test_leak_measurement.py, about 10 minutes
#   make test_fcc_slip_sys.py    one test
#   make list                    show which tests would run
#
# PYTHON can be overridden to match the interpreter pydis is installed for:
#   make PYTHON=$CONDA_PREFIX/bin/python3
#
# Each script locates the repository from its own __file__ and writes into
# output/ under the current directory.

PYTHON ?= python3
SLOW   := test_leak_measurement.py
TESTS  := $(filter-out $(SLOW), $(wildcard test_*.py))

.PHONY: all list leak clean $(TESTS) $(SLOW)

all:
	@rc=0; failed=""; \
	for t in $(TESTS); do \
	    echo "=== $$t"; \
	    MPLBACKEND=Agg $(PYTHON) $$t || { rc=1; failed="$$failed $$t"; }; \
	done; \
	if [ $$rc -eq 0 ]; then echo "all tests passed"; \
	else echo "FAILED:$$failed"; fi; \
	exit $$rc

$(TESTS) $(SLOW):
	MPLBACKEND=Agg $(PYTHON) $@

leak: $(SLOW)

list:
	@for t in $(TESTS); do echo $$t; done
	@echo "$(SLOW)  (make leak)"

clean:
	rm -rf output
