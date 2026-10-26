# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt. make and ctest must cover the same set of tests, and
# deriving both from the same rule is what keeps them from drifting apart.
#
#   make                                     run every test, non-zero status if any failed
#   make test_collision_mode_pydis_exadis.py run one test
#   make list                                show which tests would run
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=$CONDA_PREFIX/bin/python3
#
# No environment setup is needed: each script locates the repository from its own __file__.

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
