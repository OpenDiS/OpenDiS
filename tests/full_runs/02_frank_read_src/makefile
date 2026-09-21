# Run the tests in this folder.
#
# The test list is derived with $(wildcard test_*.py) so that it matches the file(GLOB ...)
# in tests/unit_tests/CMakeLists.txt. make and ctest must cover the same set of tests, and
# deriving both from the same rule is what keeps them from drifting apart.
#
#   make                      run every test, report a non-zero status if any failed
#   make test_frank_read_src_pydis_exadis_elast.py   run one test
#   make list                 show which runs would happen (a test can appear
#                             more than once, with different arguments)
#   make frank_read_src_elast_ref   regenerate the reference configuration
#
# PYTHON can be overridden to match the interpreter CMake was configured with:
#   make PYTHON=$CONDA_PREFIX/bin/python3
#
# No environment setup is needed: each script locates the repository from its own __file__.

PYTHON ?= python3
TESTS  := $(wildcard test_*.py)

# The comparison harness compares one coarsening branch per invocation, and with
# no arguments that is coarsen_mode=1. The two branches are different algorithms
# rather than settings of one, so both are covered here, as a second invocation.
# '@' stands for the space between an argument and the next, since the loop below
# splits its list on whitespace. tests/full_runs/CMakeLists.txt registers the
# same extra run with CTest, so make and ctest still cover the same set.
EXTRA_RUNS := test_frank_read_src_pydis_exadis_elast.py@--coarsen-mode@0
RUNS       := $(TESTS) $(EXTRA_RUNS)

EXAMPLES = ../../../examples/02_frank_read_src
REF_NPZ  = output/frank_read_src_elast_ref.npz

.PHONY: all list clean frank_read_src_elast_ref $(TESTS)

# Each test is run even if an earlier one fails, so one invocation reports the
# whole picture, and the names of any that failed are repeated at the end.
# Without that summary a crash is easy to miss: a segfault prints only the
# shell's own line, which names neither the test nor what happened to it, and
# it lands after a run of PASSED lines from the part of that test which did
# work. Exit status 139 is signal 11, so those are reported as signals rather
# than as an exit code nobody reads as SIGSEGV.
all:
	@rc=0; failed=""; \
	for r in $(RUNS); do \
	    t=`echo $$r | tr '@' ' '`; \
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
	@for r in $(RUNS); do echo $$r | tr '@' ' '; done

# Regenerate the reference configuration the tests here compare against.
# Either _elast example can produce it; pydis is used, and the code that wrote
# it is recorded inside the .npz as 'source'. The example writes its output
# relative to the working directory, so the .npz lands in output/ here.
frank_read_src_elast_ref:
	$(PYTHON) $(EXAMPLES)/test_frank_read_src_pydis_elast.py \
	    --no-plot --write-ref
	@echo ""
	@echo "reference written to $(REF_NPZ)"
	@echo ""
	@echo "to install it, copy it into ref_data/ alongside this makefile:"
	@echo "    mkdir -p ref_data"
	@echo "    cp $(REF_NPZ) ref_data/"
	@echo ""
	@echo "this reference encodes the current settings of the run"
	@echo "(cutoff = 0.25*Lbox, Ec = 0.0 on the exadis side, max_step = 200);"
	@echo "regenerate it whenever those change, or the tests will report a"
	@echo "legitimate improvement as a failure"

clean:
	@echo "nothing to clean"
