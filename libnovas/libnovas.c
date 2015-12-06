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

#include "Cdist/readeph0.c"
#include "Cdist/novas.c"
#include "Cdist/nutation.c"

