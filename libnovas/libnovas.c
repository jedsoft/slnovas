#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define LIBNOVAS_SOURCE

#ifndef NOVAS_SOLSYS_MODEL
# define NOVAS_SOLSYS_MODEL 3
#endif

#include "libnovas.h"
#include "_libnovas.h"

#include "Cdist/novascon.c"

#if (NOVAS_SOLSYS_MODEL == 1)
# include "Cdist/eph_manager.c"
# include "Cdist/solsys1.c"
#endif
#if (NOVAS_SOLSYS_MODEL == 2)
# include "Cdist/solsys2.c"
#endif
#if (NOVAS_SOLSYS_MODEL == 3)
# include "Cdist/solsys3.c"
#endif

#if (NOVAS_SOLSYS_MODEL != 1)
static void function_unavailable_error (const char *name, int ssv)
{
   (void) fprintf (stderr, "This version does not support %s.  Recompile with NOVAS_SOLSYS_MODEL=%d\n", name, ssv);
}
#endif

short novas_ephem_open (char *ephem_name, double *jd_begin, double *jd_end, short *denump)
{
#if (NOVAS_SOLSYS_MODEL != 1)
   (void) ephem_name; (void) jd_begin; (void) jd_end; (void) denump;
   function_unavailable_error ("novas_ephin_open", 1);
   return 1;
#else
   if (ephem_name == NULL)
     ephem_name = getenv ("JPLEPH");
   if (ephem_name == NULL)
     {
	fprintf (stderr, "novas_ephin_open: JPLEPH environemnt variable not set.\n");
	return 1;		       /* not found */
     }
   return ephem_open (ephem_name, jd_begin, jd_end, denump);
#endif
}

short novas_ephem_close (void)
{
#if (NOVAS_SOLSYS_MODEL == 1)
   return ephem_close ();
#else
   return 0;
#endif
}

#include "Cdist/readeph0.c"
#include "Cdist/novas.c"
#include "Cdist/nutation.c"

