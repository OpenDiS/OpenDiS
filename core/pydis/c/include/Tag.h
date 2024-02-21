/***************************************************************************
 *
 *  Tag.h  Define the tag struct used to uniquely identify a node as a 
 *         combination of home domain and local index in domain
 *
 ***************************************************************************/
#ifndef _Tag_h
#define _Tag_h

#include "Typedefs.h"

struct _tag {

   int domainID ;
   int index ;

} ;

#endif
