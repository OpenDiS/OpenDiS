# Run every full run, by recursing into each case folder.
#
# Case folders are discovered rather than listed, for the same reason the per-folder makefiles
# use $(wildcard test_*.py): the set of tests must match the file(GLOB ...) in CMakeLists.txt,
# and deriving both from a rule is what keeps make and ctest from drifting apart.
#
#   make                                run every test in every case folder
#   make 02_frank_read_src             run one case folder
#   make list                           show which tests would run
#   make disl_loop_elast_ref            regenerate one stored reference; every
#                                       *_ref target of every case folder can be
#                                       invoked from here
#
# PYTHON is passed down, so a single override applies everywhere:
#   make PYTHON=$CONDA_PREFIX/bin/python3
#
# ctest is the equivalent for CI: it needs a configured build tree, this does not.

PYTHON  ?= python3
SUBDIRS := $(patsubst %/,%,$(dir $(wildcard */makefile)))

# Reference-regeneration targets are defined in the case folders, one per stored
# reference, and named *_ref. They are discovered here rather than listed, for the
# same reason the test list is, and forwarded to the folder that defines them, so
# that "make disl_loop_elast_ref" works from this level and not only inside the
# case folder.
REF_TARGETS := $(shell sed -n 's/^\([a-z0-9_]*_ref\):.*/\1/p' */makefile \
                       2>/dev/null | sort -u)

.PHONY: all list clean $(SUBDIRS) $(REF_TARGETS)

# Every folder runs even if an earlier one fails, so one invocation reports the whole picture;
# the collected status is what makes this usable as a gate.
all:
	@rc=0; \
	for d in $(SUBDIRS); do \
	    echo "########## $$d"; \
	    $(MAKE) --no-print-directory -C $$d PYTHON=$(PYTHON) || rc=1; \
	done; \
	echo "##########"; \
	if [ $$rc -eq 0 ]; then echo "all full runs passed"; else echo "SOME FULL RUNS FAILED"; fi; \
	exit $$rc

$(SUBDIRS):
	@$(MAKE) --no-print-directory -C $@ PYTHON=$(PYTHON)

$(REF_TARGETS):
	@d=$$(grep -l '^$@:' */makefile | head -1 | xargs dirname); \
	echo "########## $$d"; \
	$(MAKE) --no-print-directory -C $$d PYTHON=$(PYTHON) $@

list:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$d list; done

clean:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -C $$d clean; done
