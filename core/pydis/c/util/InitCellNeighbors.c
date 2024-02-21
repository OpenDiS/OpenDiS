/****************************************************************************
 *
 *      Function:    InitCellNeighbors
 *      Description: For each native cell in domain, find all its neighbor
 *                   cells. If a neighbor cell has already been encountered
 *                   and allocated, just add its cell index to the native
 *                   cell's list of neighbors. If this is the first time,
 *                   allocate it and add it to home->cellList and
 *                   home->cellKeys, as well as to the native node's
 *                   neighbor list.
 *
 *                   If periodic boundaries are turned on and a neighbor
 *                   cell lies outside the base problem space, a cell is
 *                   allocated for the periodic image cell, but the cell
 *                   only contains offsets to apply to the corresponding
 *                   base cell, and a pointer to the base cell, which is
 *                   also allocated even if it isn't used by domain except
 *                   in shifted form. Only the base cell is added to
 *                   home->cellList
 *
 ****************************************************************************/
#include "Init.h"
#include "Home.h"
#include "Cell.h"
#include "Util.h"

void InitCellNeighbors(Home_t *home)
{
        int     iCell, ntvIdx, iX, jY, kZ, iNbr, jNbr, kNbr, nbrIdx, extend;
        int     iBase, jBase, kBase, baseIdx;
        int     nXcells, nYcells, nZcells;
        real8   xShift, yShift, zShift;
        real8   sideX, sideY, sideZ;
        Cell_t  *ntvCell, *nbrCell, *baseCell;
        Param_t *param;

        param = home->param;
   
        nXcells = param->nXcells;
        nYcells = param->nYcells;
        nZcells = param->nZcells;

        sideX = (param->maxSideX - param->minSideX);
        sideY = (param->maxSideY - param->minSideY);
        sideZ = (param->maxSideZ - param->minSideZ);

/*
 *      Loop over all cells native to the domain
 */
        for (iCell = 0; iCell < home->nativeCellCount; iCell++) {

            ntvIdx = home->cellList[iCell];
            ntvCell = home->cellKeys[ntvIdx];

/*
 *          Allocate array to hold list of neighbors of this cell
 */
            ntvCell->nbrList = (int *) malloc(26 * sizeof(int));
            ntvCell->nbrCount = 0;

/*
 *          Include all cells surrounding the native cell.  Base cells
 *          are in the range [1,nXcells] X [1,nYcells] X [1,nZcells].
 *          Cells outside this range (i.e. with any cell index 0 or
 *          n[X|Y|Z]cells) are "periodic" cells.
 */
            iX = ntvCell->xIndex;
            jY = ntvCell->yIndex;
            kZ = ntvCell->zIndex;

            for (iNbr = iX-1; iNbr <= iX+1; iNbr++) {
                for (jNbr = jY-1; jNbr <= jY+1; jNbr++) {
                    for (kNbr = kZ-1; kNbr <= kZ+1; kNbr++) {

/*
 *                      Skip the current native cell itself
 */
                        if (iNbr == iX && jNbr == jY && kNbr == kZ) continue;

/*
 *                      If neighbor cell already allocated, just add its
 *                      index to the native cell's neighbor list
 */
                        nbrIdx = EncodeCellIdx(home, iNbr, jNbr, kNbr);

                        if (home->cellKeys[nbrIdx] != 0) {
                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;
                            continue;
                        }

/*
 *                      In the following, check for cell indices outside
 *                      the base cell region. If so, and if problem is
 *                      periodic in that dimension, determine how much
 *                      the base cell must be shifted in the periodic
 *                      image, what the corresponding base cell index
 *                      is, and flag the cell as periodic. If problem
 *                      not periodic in that dimension, just skip - there
 *                      is no neighbor in that direction.
 */
                        extend = 0;        /* assume not periodic */

                        xShift = 0.0;
                        yShift = 0.0;
                        zShift = 0.0;

                        iBase = iNbr;
                        jBase = jNbr;
                        kBase = kNbr;

                        if (iNbr == 0) {    /* periodic cell */
                            if (param->xBoundType == Periodic) {
    
                                /* i-index of corresponding base cell */
                                iBase = nXcells;

                                /* correction to base X coordinate */
                                xShift = -sideX;

                                /* flag cell as "periodic" */
                                extend = 1;

                            } else continue; /* free boundary; no neighbor */
                        }

                        if (iNbr > nXcells) {    /* periodic cell */
                            if (param->xBoundType == Periodic) {

                                iBase = 1;
                                xShift = sideX;
                                extend = 1;

                            } else continue;
                        }

                        if (jNbr == 0) {    /* periodic cell */
                            if (param->yBoundType == Periodic) {

                                /* j-index of corresponding base cell */
                                jBase = nYcells;

                                /* correction to base Y coordinate */
                                yShift = -sideY;

                                /* flag cell as "periodic" */
                                extend = 1;

                            } else continue; /* free boundary; no neighbor */
                        }

                        if (jNbr > nYcells) {    /* periodic cell */
                            if (param->yBoundType == Periodic) {

                                jBase = 1;
                                yShift = sideY;
                                extend = 1;

                            } else continue;
                        }

                        if (kNbr == 0) {    /* periodic cell */
                            if (param->zBoundType == Periodic) {

                                /* k-index of corres' base cell */
                                kBase = nZcells;

                                /* correction to base Z coordinate */
                                zShift = -sideZ;

                                /* flag cell as "periodic" */
                                extend = 1;

                            } else continue; /* free boundary; no neighbor */
                        }

                        if (kNbr > nZcells) {    /* periodic cell */
                            if (param->zBoundType == Periodic) {

                                kBase = 1;
                                zShift = sideZ;
                                extend = 1;

                            } else continue;
                        }

/*
 *                      If current neighbor cell is a periodic cell,
 *                      see if base period cell has been alloc'd. If 
 *                      not, alloc now. Then alloc periodic cell and
 *                      set its shifts and a pointer to its base cell.
 *                      If neighbor is not a periodic cell, just alloc
 *                      it. In any case, add the new neighbor to the
 *                      native cell's neighbor list.
 */
                        if (extend == 0) {

/*
 *                          Not a periodic cell. Just allocate it, add
 *                          it to the domain cell list and the neighbor
 *                          list for this native cell, and put its address
 *                          in cellKeys
 */
                            nbrCell = (Cell_t *) calloc(1, sizeof(Cell_t));

                            nbrCell->xShift = xShift;
                            nbrCell->yShift = yShift;
                            nbrCell->zShift = zShift;

                            nbrCell->nbrList = 0;
                            nbrCell->baseIdx = -1;
   
/*
 *                          For base cells we need to save the per-
 *                          dimension cell indices so we don't have to make
 *                          repeated calls to DecodeCellIdx() later to figure
 *                          them out.
 */
                            nbrCell->xIndex = iNbr;
                            nbrCell->yIndex = jNbr;
                            nbrCell->zIndex = kNbr;

                            home->cellList[home->cellCount++] = nbrIdx;
                            home->cellKeys[nbrIdx] = nbrCell;

                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;

                        } else {    /*  periodic cell  */

                            baseIdx = EncodeCellIdx(home, iBase, jBase, kBase);

/*
 *                          If base cell for current periodic neighbor
 *                          not allocated, alloc now
 */
                            if (home->cellKeys[baseIdx] == 0) {

                                baseCell = (Cell_t *) calloc(1, sizeof(Cell_t));

                                baseCell->xShift = 0.0;
                                baseCell->yShift = 0.0;
                                baseCell->zShift = 0.0;

                                baseCell->nbrList = 0;
                                baseCell->baseIdx = -1;

/*
 *                              For base cells we need to save the per-
 *                              dimension cell indices so we don't have to make
 *                              repeated calls to DecodeCellIdx() later to
 *                              figure them out.
 */
                                baseCell->xIndex = iBase;
                                baseCell->yIndex = jBase;
                                baseCell->zIndex = kBase;

                                home->cellKeys[baseIdx] = baseCell;
                                home->cellList[home->cellCount++] = baseIdx;

                            }

/*
 *                          Allocate the periodic neighbor cell. NOTE
 *                          that periodic neighbor cells are not added
 *                          to the domain cell list
 */
                            nbrCell = (Cell_t *) calloc(1, sizeof(Cell_t));

                            nbrCell->xShift = xShift;
                            nbrCell->yShift = yShift;
                            nbrCell->zShift = zShift;

                            nbrCell->nbrList = 0;
                            nbrCell->domains = 0;

                            nbrCell->baseIdx = baseIdx;
/*
 *                          For ghost cells we *shouldn't* need to look up the
 *                          per-dimension cell indices, so initialize them to
 *                          invalid indices.
 */
                            nbrCell->xIndex = -1;
                            nbrCell->yIndex = -1;
                            nbrCell->zIndex = -1;

                            home->cellKeys[nbrIdx] = nbrCell;
                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;

                        }  /* end if (extend == 0)  */

                    }  /* end for (kNbr ...) */

                }  /* end for (jNbr ...)  */

            }  /* end for (iNbr ...) */

        }  /* end for (iCell ...) */

        return;
}
