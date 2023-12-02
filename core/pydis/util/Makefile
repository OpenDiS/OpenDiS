LIB_PYDIS_UTIL = pydis_util.o

all: $(LIB_PYDIS_UTIL)

Util_subset.o: Util_subset.c
	gcc -c -O3 $^ -I ../include

Stub.o: Stub.c
	gcc -c -O3 $^ -I ../include

$(LIB_PYDIS_UTIL): Util_subset.o Stub.o
	$(info --------------------------------------------------------------------)
	$(info check Stub.c for functions still need to be implemented)
	$(info --------------------------------------------------------------------)
	$(info )
	ar r $@ $^ 

clean:
	rm -f *.o
