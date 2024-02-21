/****************************************************************************
 *
 *	MirrorDomain.h:	Define the struct that holds all nodal position data
 *			from a remote (mirror) domain.  This data is used
 *			for X-window plotting and other output. Only domain 0
 *			populates and uses this structure.
 *
 ****************************************************************************/

#ifndef _MirrorDomain_h
#define _MirrorDomain_h

#include "Typedefs.h"
#include "Home.h"
#include "Tag.h"

struct _mirrordomain {
        Node_t **nodeKeys; /* indexed by node's tag.index, points */
                           /* to Node_t */
        int newNodeKeyPtr; /* elements of nodeKeys array */

        real8 *armX; /* arm[XYZ] are pointers to the corresponding */
        real8 *armY; /* X,Y and Z coordinates of nodes at the far */
        real8 *armZ; /* end of arms.  These arrays are only used */
                     /* during the process of uploading the problem*/
                     /* to task 0 when generating output, and are */
                     /* needed because some of the output functions*/
                     /* require the neighbors coordinates, but the */
                     /* entire problem space may not fit into */
                     /* task 0's memory.  So, the remote domains */
                     /* send a node's neighbor node's coords along*/
                     /* with the arms data. */
};

void FreeMirrorDomain(Home_t *home, int domIndex);

#endif
