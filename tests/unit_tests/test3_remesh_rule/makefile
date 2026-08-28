# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt, which is what keeps make and ctest covering the same
# set.
#
# Each test is run once per coarsening algorithm, since the two are different algorithms
# rather than settings of one and pydis implements both. The script's own default is
# COARSEN_MODES' last entry, 1, matching exadis and ParaDiS. CMakeLists.txt registers the
# non-default mode as its own ctest case, so the two runners still cover the same set.
#
#   make                run every test in every coarsening mode
#   make list           show which runs would happen
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=/opt/anaconda3/envs/disl_dyn/bin/python3

PYTHON ?= python3
TESTS  := $(wildcard test_*.py)
COARSEN_MODES := 0 1

.PHONY: all list clean $(TESTS)

all:
	@rc=0; \
	for t in $(TESTS); do \
	    for m in $(COARSEN_MODES); do \
	        echo "=== $$t --coarsen-mode $$m"; \
	        $(PYTHON) $$t --coarsen-mode $$m || rc=1; \
	    done; \
	done; \
	exit $$rc

$(TESTS):
	@$(PYTHON) $@

list:
	@for t in $(TESTS); do \
	    for m in $(COARSEN_MODES); do echo "$$t --coarsen-mode $$m"; done; \
	done

clean:
	@echo "nothing to clean"
