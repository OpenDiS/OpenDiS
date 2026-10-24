/* Compiles pydis' portable log/atan into the ExaDiS targets.
 *
 * Includes the implementation rather than copying it, so the two codes cannot
 * end up rounding the same function differently. It is C, and valid C++; going
 * through a .cpp means the ExaDiS targets need no C language enabled, which the
 * top-level project (LANGUAGES CXX) does not have.
 *
 * Compiled as C++, so the symbols carry C++ linkage and the declarations in
 * force_common.h alongside this file must not say extern "C".
 */
#include "../../core/pydis/c/calforce/portable_math.c"
