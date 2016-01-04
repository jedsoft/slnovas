#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define LIBNOVAS_SOURCE

#ifndef NOVAS_DATADIR
# define NOVAS_DATADIR "/usr/local/share/libnovas"
#endif

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

#ifdef fopen
# undef fopen
#endif

static char *create_datadir_filename (const char *file)
{
#ifdef NOVAS_DATADIR
   char *newfile;
   size_t len;

   len = strlen (NOVAS_DATADIR);
   newfile = (char *) malloc (len + 1 + strlen(file) + 1);
   if (newfile != NULL)
     {
	strcpy (newfile, NOVAS_DATADIR);
	newfile[len] = '/';
	strcpy (newfile + len + 1, file);
     }
   return newfile;
#else
   return NULL;
#endif
}

static FILE *_novas_fopen (const char *file, const char *mode)
{
   char *newfile;
   FILE *fp;

   if (NULL != (fp = fopen (file, mode)))
     return fp;

   if (NULL != strchr (file, '/'))
     return NULL;

   if (NULL == (newfile = create_datadir_filename (file)))
     return NULL;

   fp = fopen (newfile, mode);
   free (newfile);
   return fp;
}

