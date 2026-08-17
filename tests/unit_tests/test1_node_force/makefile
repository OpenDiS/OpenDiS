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
segseg_force_ref:
	$(PYTHON) make_parallel_ref_data.py

clean:
	@echo "nothing to clean"
