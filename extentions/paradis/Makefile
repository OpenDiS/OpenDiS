.PHONY: check-env install clean

all: check-env info install

check-env:
ifndef PARADIS_DIR
	$(error PARADIS_DIR is undefined)
endif

info:
	$(info to build the ParaDiS library)
	$(info cd $(PARADIS_DIR)/src; make libparadis.so)
	$(info cd $(PARADIS_DIR)/python; make )

install: lib/libparadis.so lib/Home.py python/paradis_util.py

lib/libparadis.so:
	ln -s $(PARADIS_DIR)/lib/libparadis.so lib

lib/Home.py:
	ln -s $(PARADIS_DIR)/lib/Home.py lib

python/paradis_util.py:
	ln -s $(PARADIS_DIR)/python/paradis_util.py python

clean:
	rm -f lib/* python/*
