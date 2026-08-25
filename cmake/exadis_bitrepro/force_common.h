/* Shim that redirects log() and atan() inside ExaDiS'
 * src/force_types/force_common.h to the portable pydis_log()/pydis_atan(),
 * without editing the ExaDiS submodule.
 *
 * WHY
 *
 * force_common.h holds SegSegForceIsotropic and its two helpers, and its eight
 * log/atan calls are the *only* thing separating ExaDiS' segment-segment force
 * from pydis' compiled ParaDiS kernel: measured over 160 pairs, the two agree
 * bit for bit once they share a libm (see
 * tests/unit_tests/test1_node_force/test_segseg_force_pydis_exadis.py). Give
 * ExaDiS the same log/atan and the two agree exactly, which is what makes a
 * reference blessed on one platform usable on another.
 *
 * HOW IT IS REACHED
 *
 * ExaDiS' force.h includes it unqualified, #include "force_common.h", and
 * force.h does not sit in the same directory, so resolution falls through to
 * the include path. Putting this directory ahead of src/force_types on that
 * path makes this file win, and #include_next picks up the real one. The
 * top-level CMakeLists.txt does that, only when PYDIS_BITREPRO_MATH is set,
 * which is only for SYS names ending in _repro.
 *
 * WHY NOT -include ON THE WHOLE TARGET
 *
 * Because "#define log" would then be live from line 1 of every translation
 * unit, and Kokkos has unqualified log() of its own in Kokkos_Complex.hpp and
 * Kokkos_MathematicalSpecialFunctions.hpp. Shadowing keeps the macros live
 * across one header instead, and they are undefined again at the bottom.
 *
 * HOW IT CAN FAIL, AND WHAT CATCHES IT
 *
 * Silently, if upstream renames the file, moves those calls elsewhere, or adds
 * an include ahead of vec.h. Nothing here would complain. What complains is the
 * test above: with ExaDiS held to zero tolerance it reads
 * "max error 0.0000e+00 ... [BITWISE]", and a detached shim turns that red. So
 * the tolerance in segseg_tables.py is part of this mechanism, not decoration.
 */
#ifndef EXADIS_BITREPRO_FORCE_COMMON_SHIM
#define EXADIS_BITREPRO_FORCE_COMMON_SHIM

/* vec.h first, so that the nested include below finds its guard already set
 * and does not get compiled under the macros. It is force_common.h's only
 * include, which is what keeps the blast radius to that one file. */
#include "vec.h"

/* No extern "C": portable_math_cxx.cpp compiles the implementation as C++,
 * so these carry C++ linkage. */
double pydis_log(double x);
double pydis_atan(double x);

#define log  pydis_log
#define atan pydis_atan

#include_next "force_common.h"

#undef log
#undef atan

#endif
