LIB_PYDIS_CALFORCE = pydis_calforce.o

# Mirrors PYDIS_BITREPRO_MATH in ../CMakeLists.txt (see that file and
# portable_math_shim.h): unset/0 by default, matching every SYS not ending
# in _repro (this Makefile has no SYS of its own; set it by hand to mirror
# one). `make PYDIS_BITREPRO_MATH=1` redirects SegSegForce.c's log()/atan()
# to the vendored, cross-machine-reproducible pydis_log()/pydis_atan()
# instead of glibc's.
PYDIS_BITREPRO_MATH ?= 0
ifeq ($(PYDIS_BITREPRO_MATH),1)
SEGSEGFORCE_FLAGS = -include portable_math_shim.h
else
SEGSEGFORCE_FLAGS =
endif

all: $(LIB_PYDIS_CALFORCE)

portable_math.o: portable_math.c
	gcc -c -O3 $^

SegSegForce.o: SegSegForce.c portable_math.h portable_math_shim.h
	gcc -c -O3 $(SEGSEGFORCE_FLAGS) $<

SegSegForce_SBN1.o: SegSegForce_SBN1.c
	gcc -c -O3 $^

SegSegForce_SBN1_SBA.o: SegSegForce_SBN1_SBA.c
	gcc -c -O3 $^

SegmentStress.o: SegmentStress.c
	gcc -c -O3 $^

StressDueToSeg.o: StressDueToSeg.c
	gcc -c -O3 $^

$(LIB_PYDIS_CALFORCE): SegSegForce.o SegSegForce_SBN1.o SegSegForce_SBN1_SBA.o SegmentStress.o StressDueToSeg.o portable_math.o
	ld -r $^ -o $@

clean:
	rm -f *.o
