SET(SOURCES 
  SegSegForce.c
  SegSegForce_SBN1.c
  SegSegForce_SBN1_SBA.c
  SegmentStress.c
  StressDueToSeg.c
)

target_sources(pydis PRIVATE ${SOURCES})

