LIB_PYDIS_COLLISION = pydis_collision.o

all: $(LIB_PYDIS_COLLISION)

GetMinDist2.o: GetMinDist2.c
	gcc -c -O3 $^

$(LIB_PYDIS_COLLISION): GetMinDist2.o
	ld -r $^ -o $@

clean:
	rm -f *.o
