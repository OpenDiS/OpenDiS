/****************************************************************************
 *
 *      Function:     InitCellNatives
 *      Description:  Find all the cells that intersect this domain.
 *                    Allocate a cell struct for each, put the address
 *                    of the cell struct at that cell's decoded index
 *                    in the cellKeys array, and add the decoded index
 *                    to the domain's cellList
 *
 ****************************************************************************/

#include <math.h>
#include "Init.h"
#include "Home.h"
#include "Cell.h"
#include "Util.h"


void InitCellNatives(Home_t *home)
{
        int     numCells, maxExtendedCells, cIndex;
        int     xIndex, yIndex, zIndex;
        int     xCellMin, yCellMin, zCellMin;
        int     xCellMax, yCellMax, zCellMax;
        real8   xMin, yMin, zMin;
        real8   xMax, yMax, zMax;
        real8   xCellSize, yCellSize, zCellSize;
        real8   cellMax, cellMin;
        Cell_t  *pCell;
        Param_t *param;

        param = home->param;

/*
 *      Determine cell dimensions
 */
        xCellSize = (param->maxSideX - param->minSideX) / param->nXcells;
        yCellSize = (param->maxSideY - param->minSideY) / param->nYcells;
        zCellSize = (param->maxSideZ - param->minSideZ) / param->nZcells;

/*
 *      Make a temporary copy of the local domain boundaries.
 */
        xMin = home->domXmin;
        xMax = home->domXmax;

        yMin = home->domYmin;
        yMax = home->domYmax;

        zMin = home->domZmin;
        zMax = home->domZmax;

/*
 *      Find the min and max cell indices (in each dimension) of cells
 *      intersecting this domain
 *
 *      FIX ME!  We should be able to replace the loops below with
 *      direct calculations of the cell indices.  Look into it...
 */
        xCellMin = 0;
        yCellMin = 0;
        zCellMin = 0;

        xCellMax = param->nXcells-1;
        yCellMax = param->nYcells-1;
        zCellMax = param->nZcells-1;

/*
 *      Get min/max cell indices in X
 */
        for (cIndex = 0; cIndex <= param->nXcells; cIndex++) {

            cellMax = rint((param->minSideX) + (cIndex * xCellSize));
            cellMin = rint((param->minSideX) + (cIndex-1) * xCellSize);

            if ((xMin >= cellMin) && (xMin < cellMax)) xCellMin = cIndex;
            if ((xMax <= cellMax) && (xMax > cellMin)) xCellMax = cIndex;
        }

/*
 *      Get min/max cell indices in Y
 */
        for (cIndex = 0; cIndex <= param->nYcells; cIndex++) {
    
            cellMax = rint((param->minSideY) + (cIndex * yCellSize));
            cellMin = rint((param->minSideY) + (cIndex-1) * yCellSize);

            if ((yMin >= cellMin) && (yMin < cellMax)) yCellMin = cIndex;
            if ((yMax <= cellMax) && (yMax > cellMin)) yCellMax = cIndex;
    }

/*
 *      Get min/max cell indices in Z
 */
        for (cIndex = 0; cIndex <= param->nZcells; cIndex++) {

            cellMax = rint((param->minSideZ) + (cIndex * zCellSize));
            cellMin = rint((param->minSideZ) + (cIndex-1) * zCellSize);

            if ((zMin >= cellMin) && (zMin < cellMax)) zCellMin = cIndex;
            if ((zMax <= cellMax) && (zMax > cellMin)) zCellMax = cIndex;
    }


/*
 *      Save the minumum and maximum cells in each dimension in param, for
 *      later use in SortNativeNodes. Save the "natural" (i.e. pre-shifted)
 *      indices. They are easier to work with in [i,j,k] space.  Shift cell
 *      indices by one to make room for periodic images at index 0
 */
        param->iCellNatMin = xCellMin-1;
        param->iCellNatMax = xCellMax-1;

        param->jCellNatMin = yCellMin-1;
        param->jCellNatMax = yCellMax-1;

        param->kCellNatMin = zCellMin-1;
        param->kCellNatMax = zCellMax-1;


/*
 *      Allocate and initialize the tag array
 */
        numCells = (param->nXcells+2) *
                   (param->nYcells+2) *
                   (param->nZcells+2);

        home->cellKeys = (Cell_t **)calloc(1, numCells * sizeof(Cell_t *));

/*
 *      Calculate the max number of cells either native to (i.e. intersecting) 
 *      this domain, or neighbor to a native cell.
 */
        maxExtendedCells = (xCellMax-xCellMin+3) *
                           (yCellMax-yCellMin+3) *
                           (zCellMax-zCellMin+3);

/*
 *      The cellList is a list of the encoded indices of each native and
 *      neighbor cell for this domain. The encoded index of the cell
 *      (X, Y, Z) is Z + (nZCells * Y) + (nZcells * nYcells * X).
 *      The native cells are always listed first in cellList; neighbor
 *      cells will be added later.
 */
        home->cellList = (int *) malloc (maxExtendedCells * sizeof(int));
        home->cellCount = 0;

/*
 *      Allocate a Cell_t struct for each native cell, and put its address in
 *      the tag array at the encoded index, and its encoded index in cellList.
 *      Since native cells cannot be periodic images, initialize periodic
 *      shifts to 0 and index of base cell to -1 (i.e. none)
 */
        for (xIndex = xCellMin; xIndex <= xCellMax; xIndex++) {
            for (yIndex = yCellMin; yIndex <= yCellMax; yIndex++) {
                for (zIndex = zCellMin; zIndex <= zCellMax; zIndex++) {
                    cIndex = EncodeCellIdx (home, xIndex, yIndex, zIndex);
                    home->cellList[home->cellCount++] = cIndex;
                    pCell = (Cell_t *) calloc (1, sizeof(Cell_t));
                    pCell->xShift = 0.0;
                    pCell->yShift = 0.0;
                    pCell->zShift = 0.0;
                    pCell->baseIdx = -1;
                    home->cellKeys[cIndex] = pCell;
/*
 *                  Also save the per-dimension indices of the cell so that we
 *                  don't have to repeatedly call DecodeCellIdx() later on.
 */
                    pCell->xIndex = xIndex;
                    pCell->yIndex = yIndex;
                    pCell->zIndex = zIndex;
                }
            }
        }

        home->nativeCellCount = home->cellCount;

/*
 *      Last thing to be done is (re)initialize some of the
 *      FM cell layers and associated info.
 */
        if (param->fmEnabled) {
            FMInit(home);
        }

        return;
}
