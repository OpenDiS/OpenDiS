#ifndef _DEBUGFUNCTIONS_H
#define _DEBUGFUNCTIONS_H
/****************************************************************************
 *
 *      DebugFunctions.h  Contains prototypes for various debug-only functions
 *
 ***************************************************************************/

void CheckSegLengths(Home_t *home, char *msg);
void CheckForNANS(Home_t *home);
void CheckForEmptySimulation(Home_t *home);
void CheckForUndefinedPlanes(Home_t *home, char *msg);

#endif
