# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt, which is what keeps make and ctest covering the same
# set.
#
#   make                run every test, report a non-zero status if any failed
#   make list           show which tests would run
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=/opt/anaconda3/bin/python3

PYTHON ?= python3
TESTS  := $(wildcard test_*.py)

.PHONY: all list clean $(TESTS)

all:
	@rc=0; \
	for t in $(TESTS); do \
	    echo "=== $$t"; \
	    $(PYTHON) $$t || rc=1; \
	done; \
	exit $$rc

$(TESTS):
	@$(PYTHON) $@

list:
	@for t in $(TESTS); do echo $$t; done

clean:
	@echo "nothing to clean"
