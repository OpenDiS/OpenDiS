LIB_PYDIS_SO = ../../lib/libpydis.so

all: $(LIB_PYDIS_SO)

util/pydis_util.o:
	cd util; make

remesh/pydis_remesh.o:
	cd remesh; make

calforce/pydis_calforce.o:
	cd calforce; make

$(LIB_PYDIS_SO): util/pydis_util.o remesh/pydis_remesh.o calforce/pydis_calforce.o
	gcc -shared $^ -o $@

clean:
	cd util; make clean
	cd remesh; make clean
	cd calforce; make clean
