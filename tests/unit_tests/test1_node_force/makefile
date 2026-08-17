# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt. make and ctest must cover the same set of tests, and
# deriving both from the same rule is what keeps them from drifting apart.
#
#   make                      run every test, report a non-zero status if any failed
#   make test_seg_stress.py   run one test
#   make list                 show which tests would run
#   make loop_node_force_ref  regenerate the nodal-force reference
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=/opt/anaconda3/bin/python3
#
# No environment setup is needed: each script locates the repository from its own __file__.

PYTHON ?= python3
TESTS  := $(wildcard test_*.py)

.PHONY: all list clean loop_node_force_ref $(TESTS)

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

# Regenerate the nodal-force reference both test_disnet_loop_force_pydis.py and
# test_disnet_loop_force_exadis.py compare against. pydis writes it and exadis only
# reads it, so there is one set of numbers rather than one per code. The script
# writes relative to the working directory, so the .npz lands in output/ here.
loop_node_force_ref:
	$(PYTHON) test_disnet_loop_force_pydis.py --write-ref

clean:
	@echo "nothing to clean"
