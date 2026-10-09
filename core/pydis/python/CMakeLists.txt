cmake_minimum_required(VERSION 3.14)

set(CALFORCE_HEADER_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../c/calforce)
# portable_math.h is here so that pydis_log/pydis_atan reach python. The
# elasticity self force has a closed form with a log in it, and matching ExaDiS'
# value of that log is what lets the two agree bitwise; numpy's log would not.
# SelfForce.h exposes ParaDiS' SelfForce/SelfForceIsotropic. portable_math.h
# exposes pydis_log/pydis_atan, which the numpy fallback for that self force
# uses so that an uncompiled run and a compiled one round the log the same way.
set(CALFORCE_HEADER_FILES SegSegForce.h SegmentStress.h StressDueToSeg.h SegSegForce_SBN1.h SegSegForce_SBN1_SBA.h SelfForce.h portable_math.h)
list(TRANSFORM CALFORCE_HEADER_FILES PREPEND ${CALFORCE_HEADER_PATH}/)

set(INCLUDE_HEADER_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../c/include)
set(INCLUDE_HEADER_FILES Home.h Init.h ParadisProto.h Util.h Force.h Timer.h OpList.h)
list(TRANSFORM INCLUDE_HEADER_FILES PREPEND ${INCLUDE_HEADER_PATH}/)

set(PYDIS_HEADERS ${CALFORCE_HEADER_FILES} ${INCLUDE_HEADER_FILES})

set(PYDIS_OPTIONS "")
set(PYDIS_OPTIONS_TO_CTYPESGEN "${PYDIS_OPTIONS}")
list(TRANSFORM PYDIS_OPTIONS_TO_CTYPESGEN PREPEND "-D;")

set(PYDIS_LIB_PY pydis_lib.py)
set(CTYPESGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ctypesgen)

# Use Clang preprocessor on macOS, fallback to GCC on other systems
if(APPLE)
  set(CC_PREPROCESS "clang -E")
else()
  set(CC_PREPROCESS "gcc -E")
endif()

message("PYDIS_HEADERS = ${PYDIS_HEADERS}")
message("CC_PREPROCESS = ${CC_PREPROCESS}")

add_custom_target(${PYDIS_LIB_PY} ALL
  COMMAND "${CTYPESGEN_DIR}/run.py" "--allow-gnu-c" "--cpp=${CC_PREPROCESS}" ${PYDIS_HEADERS} ${PYDIS_OPTIONS_TO_CTYPESGEN} -l ${LIB_PYDIS_SO} -o ${PYDIS_LIB_PY} 2>&1 | grep -v -E "ERROR:|INFO:|WARNING:" || true
  VERBATIM
)
install(FILES ${PROJECT_BINARY_DIR}/python/${PYDIS_LIB_PY} DESTINATION "${CMAKE_SOURCE_DIR}/lib")