SET(SOURCES 
  InitHome.c
  ParadisInit.c
  Param.c
  Parse.c
  DisableUnneededParams.c
  InitCellDomains.c
  InitCellNatives.c
  InitCellNeighbors.c
  Timer.c
  QueueOps.c
  SortNativeNodes.c
  Util_modified.c
  Util_subset.c
  Stub.c
)

target_sources(pydis PRIVATE ${SOURCES})

