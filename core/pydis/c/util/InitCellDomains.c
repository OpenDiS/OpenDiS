/*---------------------------------------------------------------------------
 *
 *      Function:     InitCellDomains
 *      Description:  For each cell in home->cellList, determine which
 *                    domains intersect the cell and save the list in
 *                    the corresponding cell->domains list.
 *
 *--------------------------------------------------------------------------*/
#include "Home.h"
#include "Cell.h"
#include "Decomp.h"


void InitCellDomains(Home_t *home)
{
        int     i, cellID;
        Cell_t  *cell;

        for (i = 0; i < home->cellCount; i++) {

            cellID = home->cellList[i];
            cell = home->cellKeys[cellID];

            GetCellDomainList(home, cellID, &cell->domCount, &cell->domains);
        }

        return;
}
