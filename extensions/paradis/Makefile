.PHONY: all check-env info install clean

all: check-env info install

check-env:
	@test -n "$(strip $(PARADIS_DIR))" || { echo "PARADIS_DIR is undefined" >&2; exit 1; }
	@for file in lib/libparadis.so lib/Home.py python/paradis_util.py; do \
		test -f "$(PARADIS_DIR)/$$file" || { \
			echo "Missing $(PARADIS_DIR)/$$file; build the ParaDiS library and Python wrappers first (see README.md)." >&2; \
			exit 1; \
		}; \
	done

info:
	$(info to build the ParaDiS library)
	$(info cd $(PARADIS_DIR)/src; make libparadis.so)
	$(info cd $(PARADIS_DIR)/python; make )

install: check-env lib/libparadis.so lib/Home.py python/paradis_util.py

# Each file also depends on validation so direct targets and parallel builds
# cannot create links before the source files have been checked.
lib/libparadis.so lib/Home.py python/paradis_util.py: | check-env

lib/libparadis.so:
	ln -s "$(abspath $(PARADIS_DIR))/lib/libparadis.so" "$@"

lib/Home.py:
	ln -s "$(abspath $(PARADIS_DIR))/lib/Home.py" "$@"

python/paradis_util.py:
	ln -s "$(abspath $(PARADIS_DIR))/python/paradis_util.py" "$@"

clean:
	rm -f lib/* python/*
