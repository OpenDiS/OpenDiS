# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt. make and ctest must cover the same set of tests, and
# deriving both from the same rule is what keeps them from drifting apart.
#
#   make                       run every test, report a non-zero status if any failed
#   make test_loop_topol_op_pydis.py  run one test
#   make list                  show which tests would run
#   make loop_topol_op_ref     regenerate the topology reference
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=/opt/anaconda3/bin/python3
#
# No environment setup is needed: each script locates the repository from its own __file__.

PYTHON ?= python3
TESTS  := $(wildcard test_*.py)

.PHONY: all list clean loop_topol_op_ref $(TESTS)

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

# Regenerate the network test_loop_topol_op_pydis.py compares against, i.e. the result of
# the topological operations it applies. The script writes relative to the working
# directory, so the .json lands in output/ here; installing it into ref_data/ is
# deliberate and manual, so that a change in the operations is reviewed rather than
# silently adopted.
loop_topol_op_ref:
	$(PYTHON) test_loop_topol_op_pydis.py --write-ref

clean:
	@echo "nothing to clean"
