# Run every test, by recursing into each test group (unit_tests, and full_runs once it exists).
#
# Groups are discovered rather than listed, matching the way tests/unit_tests/makefile finds
# case folders and the way each case folder finds its test_*.py. The set of tests run by make
# must stay identical to the set registered with ctest, and deriving every level from a rule
# instead of a hand-written list is what keeps the two from drifting apart.
#
#   make                                run every test in every group
#   make unit_tests                     run one group
#   make list                           show which tests would run
#   make loop_topol_op_ref              regenerate one stored reference; every
#                                       *_ref target of every case folder can be
#                                       invoked from here
#
# PYTHON is passed all the way down, so a single override applies everywhere:
#   make PYTHON=$CONDA_PREFIX/bin/python3
#
# This needs no build directory. ctest is the equivalent for CI and does need one:
#   ctest --test-dir build --output-on-failure

PYTHON ?= python3
GROUPS := $(patsubst %/,%,$(dir $(wildcard */makefile)))

# Reference-regeneration targets are defined in the case folders, one per stored
# reference, and named *_ref. They are discovered rather than listed, like the test
# groups themselves, and forwarded down through the group that holds them, so any
# stored reference can be regenerated from here.
REF_TARGETS := $(shell sed -n 's/^\([a-z0-9_]*_ref\):.*/\1/p' */*/makefile \
                       2>/dev/null | sort -u)

.PHONY: all list clean $(GROUPS) $(REF_TARGETS)

# Every group runs even if an earlier one fails, so one invocation reports the whole picture.
all:
	@rc=0; \
	for g in $(GROUPS); do \
	    echo "==================== $$g"; \
	    $(MAKE) --no-print-directory -C $$g PYTHON=$(PYTHON) || rc=1; \
	done; \
	echo "===================="; \
	if [ $$rc -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi; \
	exit $$rc

$(GROUPS):
	@$(MAKE) --no-print-directory -C $@ PYTHON=$(PYTHON)

$(REF_TARGETS):
	@g=$$(grep -l '^$@:' */*/makefile | head -1 | xargs dirname \
	      | xargs dirname); \
	$(MAKE) --no-print-directory -C $$g PYTHON=$(PYTHON) $@

list:
	@for g in $(GROUPS); do $(MAKE) --no-print-directory -C $$g list; done

clean:
	@for g in $(GROUPS); do $(MAKE) --no-print-directory -C $$g clean; done
