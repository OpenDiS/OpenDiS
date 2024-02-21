/**************************************************************************
 *
 *  Cell.h   Define the Cell struct which keeps track, for each cell, of
 *           neighboring cells, domains that intersect cell, nodes currently
 *           in cell, etc.
 *
 ***************************************************************************/

#ifndef _Cell_h
#define _Cell_h

#include "Typedefs.h"
#include "Node.h"

struct _cell {

   Node_t  *nodeQ ;    /* queue head of nodes currently in this cell */
   int      nodeCount ;/* number of nodes on nodeQ */

   int     *nbrList ;  /* list of neighbor cell encoded indices */
   int      nbrCount ; /* number of neighbor cells */

   int     *domains ;  /* domains that intersect cell (encoded indices) */
   int      domCount ; /* number of intersecting domains */

   int     baseIdx ;   /* encoded index of corresp' base cell (-1 if not */
                       /* periodic)                                      */
   real8   xShift ;    /* if periodic, amount to shift corresp' base coord */
   real8   yShift ;
   real8   zShift ;

   int     xIndex, yIndex, zIndex; /* Indices of the cell in each dimension */
                                   /* Note that these indices are in the    */
                                   /* range 0 <= xIndex <= nXcells and like-*/
                                   /* wise for Y and Z to account for the   */
                                   /* ghost cells                           */
} ;

#endif
