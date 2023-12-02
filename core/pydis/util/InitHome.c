/***************************************************************************
 *
 *	Function:	InitHome
 *	Description:	Allocate the Home_t struct for this processor, and
 *			perform any needed initializations
 *
 **************************************************************************/

#include "Init.h"
#include "Home.h"

Home_t *InitHome ()
{
	Home_t *home;

/*
 *	Calloc() zeroes allocated memory, so the only initializations
 *	needed here are any variables that should NOT start at zero!
 */
	home = (Home_t *) calloc(1, sizeof(Home_t));

	return(home);
}
