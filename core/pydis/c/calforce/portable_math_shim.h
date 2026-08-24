/* Force-included ahead of SegSegForce.c's own text (see the -include flag
 * in CMakeLists.txt / Makefile) so its log()/atan() calls resolve to the
 * portable pydis_log()/pydis_atan() from portable_math.c, without editing
 * SegSegForce.c itself.
 *
 * <math.h> must be pulled in, in full, before log/atan are redefined here:
 * glibc's math.h token-pastes the function name into its vector-ABI
 * declarations, and refuses to compile (with a #warning turned into a hard
 * error) if log or atan are already macros the first time it is processed.
 * Once math.h has run to completion under its own names, redefining the
 * two names afterward is just an ordinary textual substitution over the
 * rest of the translation unit.
 */
#ifndef PYDIS_LOG_ATAN_SHIM_H
#define PYDIS_LOG_ATAN_SHIM_H

#include <math.h>
#include "portable_math.h"

#define log pydis_log
#define atan pydis_atan

#endif
