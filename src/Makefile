LIB_OPENDIS = ../lib/libopendis.so

all: $(LIB_OPENDIS)

SegSegForce.o: SegSegForce.c
	gcc -c $^

$(LIB_OPENDIS): SegSegForce.o
	gcc -shared $^ -o $@

clean:
	rm -f *.o