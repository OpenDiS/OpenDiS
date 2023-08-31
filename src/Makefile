LIB_OPENDIS = ../lib/libopendis.so

all: $(LIB_OPENDIS)

SegSegForce.o: SegSegForce.c
	gcc -c $^

SegSegForce_SBN1.o: SegSegForce_SBN1.c
	gcc -c $^

SegSegForce_SBN1_SBA.o: SegSegForce_SBN1_SBA.c
	gcc -c $^

SegmentStress.o: SegmentStress.c
	gcc -c $^

$(LIB_OPENDIS): SegSegForce.o SegSegForce_SBN1.o SegSegForce_SBN1_SBA.o SegmentStress.o
	gcc -shared $^ -o $@

clean:
	rm -f *.o