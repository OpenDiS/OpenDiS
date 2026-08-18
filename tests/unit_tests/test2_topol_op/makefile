# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt. make and ctest must cover the same set of tests, and
# deriving both from the same rule is what keeps them from drifting apart.
#
#   make                       run every test, report a non-zero status if any failed
#   make test_loop_topol_op_pydis_exadis.py  run one test
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

# Each test is run even if an earlier one fails, so one invocation reports the
# whole picture, and the names of any that failed are repeated at the end.
# Without that summary a crash is easy to miss: a segfault prints only the
# shell's own line, which names neither the test nor what happened to it, and
# it lands after a run of PASSED lines from the part of that test which did
# work. Exit status 139 is signal 11, so those are reported as signals rather
# than as an exit code nobody reads as SIGSEGV.
all:
	@rc=0; failed=""; \
	for t in $(TESTS); do \
	    echo "=== $$t"; \
	    $(PYTHON) $$t; s=$$?; \
	    if [ $$s -ne 0 ]; then \
	        rc=1; \
	        if [ $$s -gt 128 ]; then \
	            why="killed by signal $$((s-128))"; \
	        else \
	            why="exit status $$s"; \
	        fi; \
	        echo "=== $$t FAILED, $$why"; \
	        failed="$$failed $$t"; \
	    fi; \
	done; \
	if [ -n "$$failed" ]; then \
	    echo "FAILED in $(notdir $(CURDIR)):"; \
	    for t in $$failed; do echo "    $$t"; done; \
	fi; \
	exit $$rc

$(TESTS):
	@$(PYTHON) $@

list:
	@for t in $(TESTS); do echo $$t; done

# Regenerate the network test_loop_topol_op_pydis_exadis.py compares against, i.e. the result of
# the topological operations it applies. The script writes relative to the working
# directory, so the .json lands in output/ here; installing it into ref_data/ is
# deliberate and manual, so that a change in the operations is reviewed rather than
# silently adopted.
loop_topol_op_ref:
	$(PYTHON) test_loop_topol_op_pydis_exadis.py --write-ref

clean:
	@echo "nothing to clean"
