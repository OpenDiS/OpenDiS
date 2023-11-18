LIB_OPENDIS = ../lib/libopendis.so

all: $(LIB_OPENDIS)

SegSegForce.o: SegSegForce.c
	gcc -c -O3 $^

SegSegForce_SBN1.o: SegSegForce_SBN1.c
	gcc -c -O3 $^

SegSegForce_SBN1_SBA.o: SegSegForce_SBN1_SBA.c
	gcc -c -O3 $^

SegmentStress.o: SegmentStress.c
	gcc -c -O3 $^

StressDueToSeg.o: StressDueToSeg.c
	gcc -c -O3 $^

$(LIB_OPENDIS): SegSegForce.o SegSegForce_SBN1.o SegSegForce_SBN1_SBA.o SegmentStress.o StressDueToSeg.o
	gcc -shared $^ -o $@

clean:
	rm -f *.o