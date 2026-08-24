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

.PHONY: all list clean loop_node_force_ref segseg_force_ref pydis_ref $(TESTS)

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

# Both reference targets below follow the same convention: whatever a test
# reads in goes to input_data/, whatever it checks against goes to ref_data/,
# and a new reference lands in output/ for you to inspect and copy by hand.
# Each target says where it left its files rather than leaving you to look.
#
# segseg_force_ref writes two files because the segment-pair tables are two
# things. The .dat is the fixture, 18 columns of geometry defining which pairs
# the test covers, and nothing about it needs blessing, so it goes straight to
# input_data/. The .npz is the reference, 12 columns of force, and it waits in
# output/ like every other reference here. Until 2026-08-16 these were one
# 30-column .dat sitting in ref_data/, which is why that folder still has to be
# checked for a leftover copy on an old working tree.
#
# The two are blessed at different moments, so they can fall out of step: a
# regenerated fixture is in place immediately while its reference is still in
# output/. The .npz therefore carries a digest of the geometry it was computed
# from, and the test refuses to compare rather than reporting a kernel failure
# that is really a missing copy.
#
# .gitignore carries an explicit negation for tests/*/*/input_data/*.dat and
# tests/*/*/ref_data/*.dat; the blanket *.dat rule would otherwise swallow the
# fixtures and a fresh clone could not run this test.

# Regenerate the nodal-force reference test_disnet_loop_force_pydis_exadis.py
# compares against. It is written from the pydis side and both codes are checked
# against it, so there is one set of numbers rather than one per code. The script
# writes relative to the working directory, so the .npz lands in output/ here.
loop_node_force_ref:
	$(PYTHON) test_disnet_loop_force_pydis_exadis.py --write-ref
	@echo ""
	@echo "  wrote output/loop_node_force_ref.npz"
	@echo "  inspect it, then bless it:  cp output/loop_node_force_ref.npz ref_data/"

# Regenerate the parallel / near-parallel segment-pair table. The script prints
# both paths and the copy command itself, so nothing is echoed after it here.
#
# Note what this one does that pydis_ref below does not: it rebuilds the
# *fixture* as well, and it takes the forces from the compiled ParaDiS kernel
# only after checking they agree with exadis, refusing to bless any row where
# the two disagree. Use it when the set of segment pairs changes. pydis_ref
# regenerates references for fixtures that already exist.
segseg_force_ref:
	$(PYTHON) make_parallel_ref_data.py

# Regenerate every reference in this folder from pydis, for fixtures that are
# already in place. Each test writes its own, through the same code path it
# checks, so the pydis line of each comparison reads zero against a freshly
# blessed reference and the other lines read their difference from it.
# Everything lands in output/ to be inspected and copied by hand.
pydis_ref:
	$(PYTHON) test_seg_stress_pydis.py --write-ref
	@echo ""
	$(PYTHON) test_segseg_force_pydis_exadis.py --write-ref
	@echo ""
	$(PYTHON) test_disnet_loop_force_pydis_exadis.py --write-ref

clean:
	@echo "nothing to clean"
