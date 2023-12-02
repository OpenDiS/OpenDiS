LIB_PYDIS_CALFORCE = pydis_calforce.o

all: $(LIB_PYDIS_CALFORCE)

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

$(LIB_PYDIS_CALFORCE): SegSegForce.o SegSegForce_SBN1.o SegSegForce_SBN1_SBA.o SegmentStress.o StressDueToSeg.o
	ar r $@ $^

clean:
	rm -f *.o
