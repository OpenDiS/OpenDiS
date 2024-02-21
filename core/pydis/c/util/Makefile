LIB_PYDIS_UTIL = pydis_util.o

all: $(LIB_PYDIS_UTIL)

InitHome.o: InitHome.c
	gcc -c -O3 $^ -I ../include

ParadisInit.o: ParadisInit.c
	gcc -c -O3 $^ -I ../include

Param.o: Param.c
	gcc -c -O3 $^ -I ../include

Parse.o: Parse.c
	gcc -c -O3 $^ -I ../include

DisableUnneededParams.o: DisableUnneededParams.c
	gcc -c -O3 $^ -I ../include

InitCellDomains.o: InitCellDomains.c
	gcc -c -O3 $^ -I ../include

InitCellNatives.o: InitCellNatives.c
	gcc -c -O3 $^ -I ../include

InitCellNeighbors.o: InitCellNeighbors.c
	gcc -c -O3 $^ -I ../include

Timer.o: Timer.c
	gcc -c -O3 $^ -I ../include

QueueOps.o: QueueOps.c
	gcc -c -O3 $^ -I ../include

SortNativeNodes.o: SortNativeNodes.c
	gcc -c -O3 $^ -I ../include

Util_modified.o: Util_modified.c
	gcc -c -O3 $^ -I ../include

Util_subset.o: Util_subset.c
	gcc -c -O3 $^ -I ../include

Stub.o: Stub.c
	gcc -c -O3 $^ -I ../include

$(LIB_PYDIS_UTIL): InitHome.o ParadisInit.o Param.o Parse.o DisableUnneededParams.o InitCellDomains.o InitCellNatives.o InitCellNeighbors.o Timer.o QueueOps.o SortNativeNodes.o Util_modified.o Util_subset.o Stub.o
	$(info --------------------------------------------------------------------)
	$(info check Stub.c for functions still need to be implemented)
	$(info --------------------------------------------------------------------)
	$(info )
	ld -r $^ -o $@

clean:
	rm -f *.o
