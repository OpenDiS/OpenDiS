LIB_PYDIS_REMESH = pydis_remesh.o

all: $(LIB_PYDIS_REMESH)

Remesh.o: Remesh.c
	gcc -c -O3 $^ -I ../include

RemeshRule_2.o: RemeshRule_2.c
	gcc -c -O3 $^ -I ../include

$(LIB_PYDIS_REMESH): Remesh.o RemeshRule_2.o
	ar r $@ $^ 

clean:
	rm -f *.o
