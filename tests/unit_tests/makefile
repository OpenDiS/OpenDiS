# Run every unit test, by recursing into each test case folder.
#
# Case folders are discovered rather than listed, for the same reason the per-folder makefiles
# use $(wildcard test_*.py): the set of tests must match the file(GLOB ...) in CMakeLists.txt,
# and deriving both from a rule is what keeps make and ctest from drifting apart.
#
#   make                                run every test in every case folder
#   make test1_node_force               run one case folder
#   make list                           show which tests would run
#
# PYTHON is passed down, so a single override applies everywhere:
#   make PYTHON=/opt/anaconda3/bin/python3
#
# ctest is the equivalent for CI: it needs a configured build tree, this does not.

PYTHON  ?= python3
SUBDIRS := $(patsubst %/,%,$(dir $(wildcard */makefile)))

.PHONY: all list clean $(SUBDIRS)

# Every folder runs even if an earlier one fails, so one invocation reports the whole picture;
# the collected status is what makes this usable as a gate.
all:
	@rc=0; \
	for d in $(SUBDIRS); do \
	    echo "########## $$d"; \
	    $(MAKE) --no-print-directory -C $$d PYTHON=$(PYTHON) || rc=1; \
	done; \
	echo "##########"; \
	if [ $$rc -eq 0 ]; then echo "all unit tests passed"; else echo "SOME UNIT TESTS FAILED"; fi; \
	exit $$rc

$(SUBDIRS):
	@$(MAKE) --no-print-directory -C $@ PYTHON=$(PYTHON)

list:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$d list; done

clean:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$d clean; done
