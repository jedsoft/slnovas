#include <stdio.h>
#include <stdlib.h>

#include "libnovas.h"
#include "_libnovas.h"

#define ephem_open novas_ephem_open
#define ephem_close novas_ephem_close
#include "Cdist/checkout-stars-full.c"
