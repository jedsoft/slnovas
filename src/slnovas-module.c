/* -*- mode: C; mode: fold; -*- */
/*
Copyright (C) 2016 John E. Davis

This file is part of the S-Lang NOVAS Module

The S-Lang NOVAS Module is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The S-Lang NOVAS Module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
USA.
*/

#include <stdio.h>
#include <string.h>
#include <slang.h>

#include "../libnovas/libnovas.h"

#ifdef __cplusplus
extern "C"
{
#endif
SLANG_MODULE(slnovas);
#ifdef __cplusplus
}
#endif

#include "version.h"

static SLang_Name_Type *Throw_Error_Callback;

static void throw_error (const char *func, int err)
{
   if (Throw_Error_Callback == NULL)
     {
	SLang_verror (SL_RunTime_Error, "novas error: %s: %d", func, err);
	return;
     }

   if (-1 == SLang_start_arg_list ())
     return;
   (void) SLang_push_string ((char *)func);
   (void) SLang_push_int (err);
   (void) SLang_end_arg_list ();
   if (SLang_get_error ())
     return;
   (void) SLexecute_function (Throw_Error_Callback);
}

static void interp_set_throw_callback (void)
{
   SLang_Name_Type *f;
   if (NULL == (f = SLang_pop_function ()))
     return;
   SLang_free_function (Throw_Error_Callback);
   Throw_Error_Callback = f;
}

static int push_doubles (double *a, unsigned int n)
{
   unsigned int i;

   for (i = 0; i < n; i++)
     {
	if (-1 == SLang_push_double (a[i]))
	  return -1;
     }
   return 0;
}

static int pop_doubles (double *a, unsigned int n)
{
   while (n != 0)
     {
	n--;
	if (-1 == SLang_pop_double (a+n))
	  return -1;
     }
   return 0;
}

static int pop_2double_short (double *a, double *b, short *c)
{
   if ((-1 == SLang_pop_short (c))
       || (-1 == SLang_pop_double (b))
       || (-1 == SLang_pop_double (a)))
     return -1;

   return 0;
}

static int pop_shorts (short *a, unsigned int n)
{
   while (n != 0)
     {
	n--;
	if (-1 == SLang_pop_short (a+n))
	  return -1;
     }
   return 0;
}

static int pop_doubles_from_array (double *a, unsigned int n)
{
   SLang_Array_Type *at;

   if (-1 == SLang_pop_array_of_type (&at, SLANG_DOUBLE_TYPE))
     return -1;
   if (at->num_elements != n)
     {
	SLang_verror (SL_Index_Error, "Expecting array to have %u elements, found %u",
		      n, at->num_elements);
	SLang_free_array (at);
	return -1;
     }
   memcpy (a, at->data, n*sizeof(double));
   SLang_free_array (at);
   return 0;
}

static int copy_doubles_from_array (SLang_Array_Type *at, double *a, unsigned int n)
{
   if (at->data_type != SLANG_DOUBLE_TYPE)
     {
	if (-1 == SLang_push_array (at, 0))
	  return -1;
	return pop_doubles_from_array (a, n);
     }
   if (at->num_elements != n)
     {
	SLang_verror (SL_Index_Error, "Expecting array to have %u elements, found %u",
		      n, at->num_elements);
	return -1;
     }

   memcpy (a, at->data, n*sizeof(double));

   return 0;
}

static int push_doubles_as_array (double *a, SLindex_Type n)
{
   SLang_Array_Type *at;

   at = SLang_create_array1 (SLANG_DOUBLE_TYPE, 0, NULL, &n, 1, 1);
   if (at == NULL)
     return -1;

   memcpy (at->data, a, n*sizeof(double));
   return SLang_push_array (at, 1);
}

static int pop_double_2short (double *x, short *y, short *z)
{
   if ((-1 == SLang_pop_short (z))
       || (-1 == SLang_pop_short (y))
       || (-1 == SLang_pop_double (x)))
     return -1;

   return 0;
}

typedef struct
{
   Novas_cat_entry_t cat_entry;
   char *starname;
   char *catalog;
}
Interp_Novas_cat_entry_t;
#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Interp_Novas_cat_entry_t,cname,slname,type,ro)

static SLang_CStruct_Field_Type Interp_Novas_cat_entry_Struct[] =
{
   FIELD(starname, "starname", SLANG_STRING_TYPE, 0),
   FIELD(catalog, "catalog", SLANG_STRING_TYPE, 0),
   FIELD(cat_entry.starnumber, "starnumber", SLANG_LONG_TYPE,0),
   FIELD(cat_entry.ra, "ra", SLANG_DOUBLE_TYPE,0),
   FIELD(cat_entry.dec, "dec", SLANG_DOUBLE_TYPE,0),
   FIELD(cat_entry.promora, "promora", SLANG_DOUBLE_TYPE,0),
   FIELD(cat_entry.promodec, "promodec", SLANG_DOUBLE_TYPE,0),
   FIELD(cat_entry.parallax, "parallax", SLANG_DOUBLE_TYPE,0),
   FIELD(cat_entry.radialvelocity, "radialvelocity", SLANG_DOUBLE_TYPE,0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static void free_novas_cat_entry (Interp_Novas_cat_entry_t *cat_entry)
{
   SLang_free_cstruct (cat_entry, Interp_Novas_cat_entry_Struct);
}

static int pop_novas_cat_entry (Interp_Novas_cat_entry_t *cat_entryp)
{
   if (-1 == SLang_pop_cstruct (cat_entryp, Interp_Novas_cat_entry_Struct))
     return -1;
   strncpy (cat_entryp->cat_entry.starname, cat_entryp->starname, NOVAS_SIZE_OF_OBJ_NAME);
   cat_entryp->cat_entry.starname[NOVAS_SIZE_OF_OBJ_NAME-1] = 0;
   strncpy (cat_entryp->cat_entry.catalog, cat_entryp->catalog, NOVAS_SIZE_OF_CAT_NAME);
   cat_entryp->cat_entry.catalog[NOVAS_SIZE_OF_CAT_NAME] = 0;

   return 0;
}

static int push_novas_cat_entry (Interp_Novas_cat_entry_t *cat_entry)
{
   int ret;

   if (NULL == (cat_entry->starname = SLang_create_slstring (cat_entry->cat_entry.starname)))
     return -1;

   if (NULL == (cat_entry->catalog = SLang_create_slstring (cat_entry->cat_entry.catalog)))
     {
	SLfree (cat_entry->starname);
	return -1;
     }

   ret = SLang_push_cstruct (cat_entry, Interp_Novas_cat_entry_Struct);
   free_novas_cat_entry (cat_entry);
   return ret;
}

typedef struct
{
   Novas_object_t object;
   SLang_Struct_Type *star;
   char *name;
}
Interp_Novas_object_t;

#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Interp_Novas_object_t,cname,slname,type,ro)

static SLang_CStruct_Field_Type Interp_Novas_object_Struct[] =
{
   FIELD(object.type, "type", SLANG_SHORT_TYPE, 0),
   FIELD(object.number, "number", SLANG_SHORT_TYPE, 0),
   FIELD(name, "name", SLANG_STRING_TYPE,0),
   FIELD(star, "star", SLANG_STRUCT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static void free_novas_object (Interp_Novas_object_t *obj)
{
   SLang_free_cstruct (obj, Interp_Novas_object_Struct);
}

static int pop_novas_object (Interp_Novas_object_t *obj)
{
   Interp_Novas_cat_entry_t star;

   if (-1 == SLang_pop_cstruct (obj, Interp_Novas_object_Struct))
     return -1;
   strncpy (obj->object.name, obj->name, NOVAS_SIZE_OF_OBJ_NAME);
   obj->object.name[NOVAS_SIZE_OF_OBJ_NAME-1] = 0;

   if ((-1 == SLang_push_struct (obj->star))
       || (-1 == pop_novas_cat_entry (&star)))
     {
	free_novas_object (obj);
	return -1;
     }
   obj->object.star = star.cat_entry;
   free_novas_cat_entry (&star);

   return 0;
}


#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Novas_on_surface_t,cname,slname,type,ro)
static SLang_CStruct_Field_Type Interp_Novas_on_surface_Struct[] =
{
   FIELD(latitude, "latitude", SLANG_DOUBLE_TYPE, 0),
   FIELD(longitude, "longitude", SLANG_DOUBLE_TYPE, 0),
   FIELD(height, "height", SLANG_DOUBLE_TYPE, 0),
   FIELD(temperature, "temperature", SLANG_DOUBLE_TYPE, 0),
   FIELD(pressure, "pressure", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static int pop_novas_on_surface (Novas_on_surface_t *surface)
{
   return SLang_pop_cstruct (surface, Interp_Novas_on_surface_Struct);
}

typedef struct
{
   Novas_in_space_t in_space;
   SLang_Array_Type *sc_pos;
   SLang_Array_Type *sc_vel;
}
Interp_Novas_in_space_t;

#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Interp_Novas_in_space_t,cname,slname,type,ro)
static SLang_CStruct_Field_Type Interp_Novas_in_space_Struct[] =
{
   FIELD(sc_pos, "sc_pos", SLANG_ARRAY_TYPE, 0),
   FIELD(sc_vel, "sc_vel", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static int pop_novas_in_space (Interp_Novas_in_space_t *inspace)
{
   if (-1 == SLang_pop_cstruct (inspace, Interp_Novas_in_space_Struct))
     return -1;

   if ((-1 == copy_doubles_from_array (inspace->sc_pos, inspace->in_space.sc_pos, 3))
       || (-1 == copy_doubles_from_array (inspace->sc_vel, inspace->in_space.sc_vel, 3)))
     {
	SLang_free_cstruct (inspace, Interp_Novas_in_space_Struct);
	return -1;
     }

   return 0;
}

typedef struct
{
   Novas_observer_t observer;
   SLang_Struct_Type *on_surf;
   SLang_Struct_Type *near_earth;
}
Interp_Novas_observer_t;
#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Interp_Novas_observer_t,cname,slname,type,ro)
static SLang_CStruct_Field_Type Interp_Novas_observer_Struct [] =
{
   FIELD(observer.where, "where", SLANG_SHORT_TYPE, 0),
   FIELD(on_surf, "on_surf", SLANG_STRUCT_TYPE, 0),
   FIELD(near_earth, "near_earth", SLANG_STRUCT_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static void free_novas_observer (Interp_Novas_observer_t *observer)
{
   SLang_free_cstruct (observer, Interp_Novas_observer_Struct);
}

static int pop_novas_observer (Interp_Novas_observer_t *observer)
{
   Interp_Novas_in_space_t near_earth;
   Novas_on_surface_t on_surf;

   if (-1 == SLang_pop_cstruct (observer, Interp_Novas_observer_Struct))
     return -1;
   if ((-1 == SLang_push_struct (observer->on_surf))
       || (-1 == pop_novas_on_surface (&on_surf)))
     {
	free_novas_observer (observer);
	return -1;
     }
   observer->observer.on_surf = on_surf;
   SLang_free_cstruct (&on_surf, Interp_Novas_on_surface_Struct);

   if ((-1 == SLang_push_struct (observer->near_earth))
       || (-1 == pop_novas_in_space (&near_earth)))
     {
	SLang_free_cstruct (observer, Interp_Novas_observer_Struct);
	return -1;
     }
   observer->observer.near_earth = near_earth.in_space;
   SLang_free_cstruct (&near_earth, Interp_Novas_in_space_Struct);

   return 0;
}

typedef struct
{
   Novas_sky_pos_t sky_pos;
   SLang_Array_Type *r_hat;
}
Interp_Novas_sky_pos_t;

#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Interp_Novas_sky_pos_t,cname,slname,type,ro)
static SLang_CStruct_Field_Type Interp_Novas_sky_pos_Struct [] =
{
   FIELD(sky_pos.ra, "ra", SLANG_DOUBLE_TYPE, 0),
   FIELD(sky_pos.dec, "dec", SLANG_DOUBLE_TYPE, 0),
   FIELD(sky_pos.dis, "dis", SLANG_DOUBLE_TYPE, 0),
   FIELD(sky_pos.rv, "rv", SLANG_DOUBLE_TYPE, 0),
   FIELD(r_hat, "r_hat", SLANG_ARRAY_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static void free_novas_sky_pos (Interp_Novas_sky_pos_t *sky_pos)
{
   SLang_free_cstruct (sky_pos, Interp_Novas_sky_pos_Struct);
}

static int push_novas_sky_pos (Interp_Novas_sky_pos_t *sky_pos)
{
   SLindex_Type dim0 = 3;
   int ret;

   if (NULL == (sky_pos->r_hat = SLang_create_array1 (SLANG_DOUBLE_TYPE, 0, NULL, &dim0, 1, 1)))
     return -1;

   memcpy (sky_pos->r_hat->data, sky_pos->sky_pos.r_hat, dim0*sizeof(double));
   ret = SLang_push_cstruct (sky_pos, Interp_Novas_sky_pos_Struct);
   free_novas_sky_pos (sky_pos);
   return ret;
}

#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(Novas_ra_of_cio_t,cname,slname,type,ro)
static SLang_CStruct_Field_Type Interp_Novas_ra_of_cio_Struct [] =
{
   FIELD(jd_tdb, "jd_tdb", SLANG_DOUBLE_TYPE, 0),
   FIELD(ra_cio, "ra_cio", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD

static int interp_xxx_star (const char *fname,
			      short (*func)(double, Novas_cat_entry_t*, short, double*, double*))
{
   double radec[2], jd_tt;
   Interp_Novas_cat_entry_t cat_entry;
   int ret;
   short accuracy;

   if ((-1 == SLang_pop_short (&accuracy))
       || (-1 == pop_novas_cat_entry (&cat_entry)))
     return -1;

   if (-1 == SLang_pop_double (&jd_tt))
     {
	SLang_free_cstruct (&cat_entry, Interp_Novas_cat_entry_Struct);
	return -1;
     }

   ret = (*func)(jd_tt, &cat_entry.cat_entry, accuracy, radec, radec+1);
   SLang_free_cstruct (&cat_entry, Interp_Novas_cat_entry_Struct);

   if (ret == 0)
     return push_doubles (radec, 2);

   throw_error (fname, ret);
   return -1;
}

static int interp_xxx_planet (const char *fname,
			      short (*func)(double, Novas_object_t*, short, double*, double*, double *))
{
   double radecdis[3], jd_tt;
   Interp_Novas_object_t obj;
   int ret;
   short accuracy;

   if ((-1 == SLang_pop_short (&accuracy))
       || (-1 == pop_novas_object (&obj)))
     return -1;

   if (-1 == SLang_pop_double (&jd_tt))
     {
	SLang_free_cstruct (&obj, Interp_Novas_object_Struct);
	return -1;
     }

   ret = (*func)(jd_tt, &obj.object, accuracy, radecdis, radecdis+1, radecdis+2);
   SLang_free_cstruct (&obj, Interp_Novas_object_Struct);

   if (ret == 0)
     return push_doubles (radecdis, 3);

   throw_error (fname, ret);
   return -1;
}

static int interp_topolocal_star (const char *fname,
				  short (*func)(double, double, Novas_cat_entry_t*, Novas_on_surface_t*, short,
						double *, double *))
{
   double radec[2], jd_tt, delta_t;
   Interp_Novas_cat_entry_t cat_entry;
   Novas_on_surface_t surface;
   int ret;
   short accuracy;

   if ((-1 == SLang_pop_short (&accuracy))
       || (-1 == pop_novas_on_surface (&surface))
       || (-1 == pop_novas_cat_entry (&cat_entry)))
     return -1;

   if ((-1 == SLang_pop_double (&delta_t))
       || (-1 == SLang_pop_double (&jd_tt)))
     {
	free_novas_cat_entry (&cat_entry);
	return -1;
     }

   ret = (*func)(jd_tt, delta_t, &cat_entry.cat_entry, &surface, accuracy,
		 radec, radec+1);

   free_novas_cat_entry (&cat_entry);

   if (ret == 0)
     return push_doubles (radec, 2);

   throw_error (fname, ret);
   return -1;
}

static int interp_topolocal_planet (const char *fname,
				    short (*func)(double, Novas_object_t*, double,
						  Novas_on_surface_t*, short,
						  double *, double *, double *))
{
   double radecdis[3], jd_tt, delta_t;
   Interp_Novas_object_t ss_body;
   Novas_on_surface_t surface;
   int ret;
   short accuracy;

   if ((-1 == SLang_pop_short (&accuracy))
       || (-1 == pop_novas_on_surface (&surface))
       || (-1 == SLang_pop_double (&delta_t))
       || (-1 == pop_novas_object (&ss_body)))
     return -1;

   if (-1 == SLang_pop_double (&jd_tt))
     {
	free_novas_object (&ss_body);
	return -1;
     }

   ret = (*func)(jd_tt, &ss_body.object, delta_t, &surface, accuracy,
		 radecdis, radecdis+1, radecdis+2);

   free_novas_object (&ss_body);

   if (ret == 0)
     return push_doubles (radecdis, 3);

   throw_error (fname, ret);
   return -1;
}

static void interp_app_star (void)
{
   (void) interp_xxx_star ("app_star", novas_app_star);
}

static void interp_virtual_star (void)
{
   (void) interp_xxx_star ("virtual_star", novas_virtual_star);
}

static void interp_astro_star (void)
{
   (void) interp_xxx_star ("astro_star", novas_astro_star);
}

static void interp_app_planet (void)
{
   (void) interp_xxx_planet ("app_planet", novas_app_planet);
}

static void interp_virtual_planet (void)
{
   (void) interp_xxx_planet ("virtual_planet", novas_virtual_planet);
}

static void interp_astro_planet (void)
{
   (void) interp_xxx_planet ("astro_planet", novas_astro_planet);
}

static void interp_topo_star (void)
{
   (void) interp_topolocal_star ("topo_star", novas_topo_star);
}

static void interp_local_star (void)
{
   (void) interp_topolocal_star ("local_star", novas_local_star);
}

static void interp_local_planet (void)
{
   (void) interp_topolocal_planet ("local_planet", novas_local_planet);
}

static void interp_topo_planet (void)
{
   (void) interp_topolocal_planet ("topo_planet", novas_topo_planet);
}

static void interp_mean_star (void)
{
   double ttradec[3];
   double radec[2];
   int ret;
   short accuracy;

   if ((-1 == SLang_pop_short (&accuracy))
       || (-1 == pop_doubles (ttradec, 3)))
     return;

   ret = novas_mean_star (ttradec[0], ttradec[1], ttradec[2], accuracy,
			  radec, radec+1);
   if (ret != 0)
     {
	throw_error ("mean_star", ret);
	return;
     }
   (void) push_doubles (radec, 2);
}

static void interp_place (void)
{
   double jd_tt, delta_t;
   Interp_Novas_object_t cel_object;
   Interp_Novas_observer_t location;
   Interp_Novas_sky_pos_t output;
   int ret;
   short coord_sys, accuracy;

   if ((-1 == pop_double_2short (&delta_t, &coord_sys, &accuracy))
       || (-1 == pop_novas_observer (&location)))
     return;

   if (-1 == pop_novas_object (&cel_object))
     {
	free_novas_object (&cel_object);
	return;
     }
   if (0 == SLang_pop_double (&jd_tt))
     {
	ret = novas_place (jd_tt, &cel_object.object, &location.observer, delta_t,
			   coord_sys, accuracy, &output.sky_pos);
     }
   free_novas_object (&cel_object);
   free_novas_observer (&location);

   if (ret != 0)
     {
	throw_error ("place", ret);
	return;
     }

   (void) push_novas_sky_pos (&output);
}

static void interp_equ2gal (double *rai, double *deci)
{
   double xy[2];

   novas_equ2gal (*rai, *deci, xy, xy+1);
   (void) push_doubles (xy, 2);
}

static void interp_equ2ecl (double *jd_tt, short *coord_sys, short *accuracy, double *ra, double *dec)
{
   int ret;
   double elonlat[2];

   ret = novas_equ2ecl (*jd_tt, *coord_sys, *accuracy, *ra, *dec, elonlat, elonlat+1);
   if (ret == 0)
     {
	(void) push_doubles (elonlat, 2);
	return;
     }
   throw_error ("equ2ecl", ret);
}

static double interp_julian_date (void)
{
   short ymd[3];
   double hour;

   if (-1 == SLang_pop_double (&hour))
     return -1.0;
   if (-1 == pop_shorts (ymd, 3))
     return -1.0;
   return novas_julian_date (ymd[0], ymd[1], ymd[2], hour);
}

static void interp_equ2hor (void)
{
   double jd_ut1, delta_t, ra, dec;
   Novas_on_surface_t location;
   double zz[4], xpyp[2];
   short accuracy, ref_option;

   if (-1 == pop_2double_short (&ra, &dec, &ref_option))
     return;
   if (-1 == pop_novas_on_surface (&location))
     return;
   if (-1 == pop_doubles (xpyp, 2))
     return;
   if (-1 == pop_2double_short (&jd_ut1, &delta_t, &accuracy))
     return;

   novas_equ2hor (jd_ut1, delta_t, accuracy, xpyp[0], xpyp[1], &location, ra, dec, ref_option,
		  zz, zz+1, zz+2, zz+3);

   (void) push_doubles (zz, 4);
}

static void interp_sidereal_time (double *jd_high, double *jd_low, double *delta_t,
				  short *gst_type, short *method, short *accuracy)
{
   double gst;
   int ret;

   ret = novas_sidereal_time(*jd_high, *jd_low, *delta_t, *gst_type, *method, *accuracy, &gst);
   if (ret == 0)
     {
	(void) SLang_push_double (gst);
	return;
     }
   throw_error ("sidereal_time", ret);
}

static void interp_era (double *jd_high, double *jd_low)
{
   double x = novas_era (*jd_high, *jd_low);
   (void) SLang_push_double (x);
}

static void interp_ephemeris (void)
{
   Interp_Novas_object_t cel_obj;
   double jd[2], pos[3], vel[3];
   short oa[2], ret;

   if (-1 == pop_shorts (oa, 2))
     return;
   if (-1 == pop_novas_object (&cel_obj))
     return;
   if (-1 == pop_doubles_from_array (jd, 2))
     {
	free_novas_object (&cel_obj);
	return;
     }

   ret = novas_ephemeris (jd, &cel_obj.object, oa[0], oa[1], pos, vel);
   free_novas_object (&cel_obj);

   if (ret == 0)
     {
	(void) push_doubles_as_array (pos, 3);
	(void) push_doubles_as_array (vel, 3);
	return;
     }

   throw_error ("ephemeris", ret);
}

static void interp_equ2exx_vec (const char *name,
				short (*func)(double, short, short, double *, double *))
{
   double jdtt;
   double pos1[3], pos2[3];
   int ret;
   short coordsys, accuracy;

   if (-1 == pop_doubles_from_array (pos1, 3))
     return;
   if (-1 == pop_double_2short (&jdtt, &coordsys, &accuracy))
     return;

   ret = (*func) (jdtt, coordsys, accuracy, pos1, pos2);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos2, 3);
	return;
     }

   throw_error (name, ret);
}

static void interp_equ2ecl_vec (void)
{
   interp_equ2exx_vec ("equ2ecl_vec", novas_equ2ecl_vec);
}

static void interp_ecl2equ_vec (void)
{
   interp_equ2exx_vec ("equ2equ_vec", novas_ecl2equ_vec);
}

static void interp_gcrs2equ (double *jdtt, short *coordsys, short *acc, double *rag, double *decg)
{
   double radec[2];
   int ret;

   ret = novas_gcrs2equ (*jdtt, *coordsys, *acc, *rag, *decg, radec, radec+1);
   if (ret == 0)
     {
	(void) push_doubles (radec, 2);
	return;
     }
   throw_error ("gcrs2equ", ret);
}

static void interp_vector2radec (void)
{
   double pos[3];
   double radec[2];
   int ret;

   if (-1 == pop_doubles_from_array (pos, 3))
     return;

   ret = novas_vector2radec (pos, radec, radec+1);
   if (ret == 0)
     {
	(void) push_doubles (radec, 2);
	return;
     }

   throw_error ("vector2radec", ret);
}

static void interp_ter_cel (const char *fname,
			    short (*f)(double, double, double, short, short, short,
				       double, double, double*, double*))
{
   double hld[3], xpyp[2];
   double vec1[3], vec2[3];
   short mao[3];
   int ret;

   if ((-1 == pop_doubles_from_array (vec1, 3))
       || (-1 == pop_doubles (xpyp, 2))
       || (-1 == pop_shorts (mao, 3))
       || (-1 == pop_doubles (hld, 3)))
     return;

   ret = (*f) (hld[0], hld[1], hld[2], mao[0], mao[1], mao[2],
	       xpyp[0], xpyp[1], vec1, vec2);
   if (ret == 0)
     {
	(void) push_doubles_as_array (vec2, 3);
	return;
     }

   throw_error (fname, ret);
}

static void interp_cel2ter (void)
{
   interp_ter_cel ("cel2ter", novas_cel2ter);
}

static void interp_ter2cel (void)
{
   interp_ter_cel ("ter2cel", novas_ter2cel);
}


static void interp_ephem_open (char *file)
{
   double jd2[2];
   short de;

   int ret = novas_ephem_open (file, jd2, jd2+1, &de);

   if (ret == 0)
     {
	(void) push_doubles (jd2, 2);
	(void) SLang_push_short (de);
	return;
     }

   throw_error ("ephem_open", ret);
}

static void interp_ephem_close (void)
{
   novas_ephem_close ();
}

static void interp_aberration (void)
{
   double pos[3], pos2[3], ve[3];
   double lt;

   if ((-1 == SLang_pop_double (&lt))
       || (-1 == pop_doubles_from_array (ve, 3))
       || (-1 == pop_doubles_from_array (pos, 3)))
     return;

   novas_aberration (pos, ve, lt, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

/* (pos2, light_time) = bary2obs (pos, pos_obs) */
static void interp_bary2obs (void)
{
   double pos[3], pos_obs[3], pos2[3], light_time;

   if ((-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_doubles_from_array (pos, 3)))
     return;

   novas_bary2obs (pos, pos_obs, pos2, &light_time);
   (void) push_doubles_as_array (pos2, 3);
   (void) SLang_push_double (light_time);
}

typedef struct
{
   short year;
   short month;
   short day;
   double hour;
}
Cal_Date_Type;
static SLang_CStruct_Field_Type Cal_Date_Struct[] =
{
   MAKE_CSTRUCT_FIELD(Cal_Date_Type, year, "year", SLANG_SHORT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cal_Date_Type, month, "month", SLANG_SHORT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cal_Date_Type, day, "day", SLANG_SHORT_TYPE, 0),
   MAKE_CSTRUCT_FIELD(Cal_Date_Type, hour, "hour", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};

static void interp_cal_date (double *tdj)
{
   Cal_Date_Type c;

   novas_cal_date (*tdj, &c.year, &c.month, &c.day, &c.hour);

   (void) SLang_push_cstruct (&c, Cal_Date_Struct);
}

static void interp_cel_pole (double *t, short *type, double *d1, double *d2)
{
   int ret = novas_cel_pole (*t, *type, *d1, *d2);
   if (ret == 0)
     return;
   throw_error ("cel_pole", ret);
}

/* (jd_tdb[], ra_cio[]) = cio_array (jd_tdb, npts); */
static void interp_cio_array (double *jd_tdb, long *n)
{
   SLang_Array_Type *at_t, *at_ra;
   double *t, *ra;
   Novas_ra_of_cio_t *cio;
   SLindex_Type i, npts;
   int ret;

   npts = (SLindex_Type) *n;
   if ((long) npts != *n)
     {
	SLang_verror (SL_LimitExceeded_Error, "Array size larger than supported pass to cio_array");
	return;
     }
   if (NULL == (cio = (Novas_ra_of_cio_t *)SLmalloc (sizeof(Novas_ra_of_cio_t) * npts)))
     return;

   ret = novas_cio_array (*jd_tdb, npts, cio);
   if (ret != 0)
     {
	throw_error ("cio_array", ret);
	goto free_and_return;
	return;
     }
   if (NULL == (at_t = SLang_create_array1 (SLANG_DOUBLE_TYPE, 0, NULL, &npts, 1, 1)))
     goto free_and_return;
   if (NULL == (at_ra = SLang_create_array1 (SLANG_DOUBLE_TYPE, 0, NULL, &npts, 1, 1)))
     {
	SLang_free_array (at_t);
	goto free_and_return;
     }

   t = (double *)at_t->data;
   ra = (double *)at_ra->data;

   for (i = 0; i < npts; i++)
     {
	t[i] = cio[i].jd_tdb;
	ra[i] = cio[i].ra_cio;
     }

   (void) SLang_push_array (at_t, 1);
   (void) SLang_push_array (at_ra, 1);
   /* drop */
free_and_return:
   SLfree ((char *)cio);
}

static void interp_cio_basis (double *t, double *ra, short *ref, short *acc)
{
   double x[3], y[3], z[3];
   int ret;

   ret = novas_cio_basis (*t, *ra, *ref, *acc, x, y, z);
   if (ret == 0)
     {
	(void) push_doubles_as_array (x, 3);
	(void) push_doubles_as_array (y, 3);
	(void) push_doubles_as_array (z, 3);
	return;
     }
   throw_error ("cio_basis", ret);
}

/* (ra_cio, ref_sys) = cio_location (jd_tdb, accuracy); */
static void interp_cio_location (double *tdb, short *accuracy)
{
   double ra_cio;
   int ret;
   short ref_sys;

   ret = novas_cio_location (*tdb, *accuracy, &ra_cio, &ref_sys);
   if (ret == 0)
     {
	(void) SLang_push_double (ra_cio);
	(void) SLang_push_short (ref_sys);
	return;
     }
   throw_error ("cio_location", ret);
}

static double interp_cio_ra (double *tt, short *acc)
{
   double ra;
   int ret;

   ret = novas_cio_ra (*tt, *acc, &ra);
   if (ret == 0)
     return ra;
   throw_error ("cio_ra", ret);
   return -1.0;
}

static double interp_d_light (void)
{
   double pos1[3], pos_obs[3];

   if ((-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_doubles_from_array (pos1, 3)))
     return -1.0;

   return novas_d_light (pos1, pos_obs);
}

static double interp_ee_ct (double *thi, double *tlo, short *acc)
{
   return novas_ee_ct (*thi, *tlo, *acc);
}

typedef struct
{
   double mean_obliquity;
   double true_obliquity;
   double ee;
   double psi;
   double eps;
}
E_Tilt_Type;
#define FIELD(cname,slname,type,ro) \
   MAKE_CSTRUCT_FIELD(E_Tilt_Type,cname,slname,type,ro)
static SLang_CStruct_Field_Type E_Tilt_Struct[] =
{
   FIELD(true_obliquity, "mean_obliquity", SLANG_DOUBLE_TYPE, 0),
   FIELD(true_obliquity, "true_obliquity", SLANG_DOUBLE_TYPE, 0),
   FIELD(ee, "ee", SLANG_DOUBLE_TYPE, 0),
   FIELD(psi, "psi", SLANG_DOUBLE_TYPE, 0),
   FIELD(eps, "eps", SLANG_DOUBLE_TYPE, 0),
   SLANG_END_CSTRUCT_TABLE
};
#undef FIELD
static void interp_e_tilt (double *tdb, short *acc)
{
   E_Tilt_Type e;
   novas_e_tilt (*tdb, *acc, &e.mean_obliquity, &e.true_obliquity, &e.ee,
		 &e.psi, &e.eps);
   (void) SLang_push_cstruct (&e, E_Tilt_Struct);
}

static void interp_frame_tie (void)
{
   double pos1[3], pos2[3];
   short dir;

   if ((-1 == SLang_pop_short (&dir))
       || (-1 == pop_doubles_from_array (pos1, 3)))
     return;

   novas_frame_tie (pos1, dir, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

static void interp_fund_args (double *t)
{
   double a[5];
   novas_fund_args (*t, a);
   (void) push_doubles_as_array (a, 5);
}

static void interp_geo_posvel (void)
{
   double tt, dt, pos[3], vel[3];
   Interp_Novas_observer_t obs;
   int ret;
   short acc;

   if (-1 == pop_novas_observer (&obs))
     return;

   if (-1 == pop_2double_short (&tt, &dt, &acc))
     {
	free_novas_observer (&obs);
	return;
     }
   ret = novas_geo_posvel (tt, dt, acc, &obs.observer, pos, vel);
   free_novas_observer (&obs);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos, 3);
	(void) push_doubles_as_array (vel, 3);
	return;
     }
   throw_error ("geo_posvel", ret);
}

static void interp_grav_def (void)
{
   double tdb, pos1[3], pos_obs[3], pos2[3];
   int ret;
   short loc, acc;

   if ((-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_doubles_from_array (pos1, 3))
       || (-1 == pop_double_2short (&tdb, &loc, &acc)))
     return;

   ret = novas_grav_def (tdb, loc, acc, pos1, pos_obs, pos2);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos2, 3);
	return;
     }
   throw_error ("grav_def", ret);
}

static void interp_grav_vec (void)
{
   double pos1[3], pos2[3], pos_obs[3], pos_body[3], rmass;

   if ((-1 == SLang_pop_double (&rmass))
       || (-1 == pop_doubles_from_array (pos_body, 3))
       || (-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_doubles_from_array (pos1, 3)))
     return;

   novas_grav_vec (pos1, pos_obs, pos_body, rmass, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

static double interp_ira_equinox (double *tdb, short *equinox, short *acc)
{
   return novas_ira_equinox (*tdb, *equinox, *acc);
}

/* (pos, tlight) = light_time (tdb, object, pos_obs, tlight0, accuracy); */
static void interp_light_time (void)
{
   double tdb, pos_obs[3], pos[3], tlight0, tlight;
   Interp_Novas_object_t obj;
   int ret;
   short acc;

   if ((-1 == SLang_pop_short (&acc))
       || (-1 == SLang_pop_double (&tlight0))
       || (-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_novas_object (&obj)))
     return;

   if (-1 == SLang_pop_double (&tdb))
     {
	free_novas_object (&obj);
	return;
     }

   ret = novas_light_time (tdb, &obj.object, pos_obs, tlight0, acc, pos, &tlight);
   free_novas_object (&obj);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos, 3);
	(void) SLang_push_double (tlight);
	return;
     }

   throw_error ("light_time", ret);
}

static void interp_limb_angle (void)
{
   double pos_obj[3], pos_obs[3], angles[2];

   if ((-1 == pop_doubles_from_array (pos_obs, 3))
       || (-1 == pop_doubles_from_array (pos_obj, 3)))
     return;

   novas_limb_angle (pos_obj, pos_obs, angles, angles+1);
   (void) push_doubles (angles, 2);
}

static double interp_mean_obliq (double *tdb)
{
   return novas_mean_obliq (*tdb);
}

static double interp_norm_ang (double *a)
{
   return novas_norm_ang (*a);
}

static void interp_nutation (void)
{
   double tdb, pos[3], pos2[3];
   short acc, dir;

   if ((-1 == pop_doubles_from_array (pos, 3))
       || (-1 == pop_double_2short (&tdb, &dir, &acc)))
     return;
   novas_nutation (tdb, dir, acc, pos, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

/* (dpsi, deps) = nutation_angles (t, accuracy); */
static void interp_nutation_angles (double *t, short *acc)
{
   double angles[2];

   novas_nutation_angles (*t, *acc, angles, angles+1);
   (void) push_doubles (angles, 2);
}

static void interp_planet_ephemeris (void)
{
   double tjd[2], pos[3], vel[3];
   int ret;
   short tc[2];

   if ((-1 == pop_shorts (tc, 2))
       || (-1 == pop_doubles_from_array (tjd, 2)))
     return;

   ret = novas_planet_ephemeris (tjd, tc[0], tc[1], pos, vel);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos, 3);
	(void) push_doubles_as_array (vel, 3);
	return;
     }
   throw_error ("planet_ephemeris", ret);
}

static void interp_precession (void)
{
   double tdb1, tdb2, pos1[3], pos2[3];
   int ret;

   if ((-1 == SLang_pop_double (&tdb2))
       || (-1 == pop_doubles_from_array (pos1, 3))
       || (-1 == SLang_pop_double (&tdb1)))
     return;

   ret = novas_precession (tdb1, pos1, tdb2, pos2);
   if (ret == 0)
     {
	(void) push_doubles_as_array (pos2, 3);
	return;
     }

   throw_error ("precession", ret);
}

static void interp_proper_motion (void)
{
   double tdb1, tdb2, pos[3], vel[3], pos2[3];

   if ((-1 == SLang_pop_double (&tdb2))
       || (-1 == pop_doubles_from_array (vel, 3))
       || (-1 == pop_doubles_from_array (pos, 3))
       || (-1 == SLang_pop_double (&tdb1)))
     return;

   novas_proper_motion (tdb1, pos, vel, tdb2, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

static void interp_radec2vector (double *ra, double *dec, double *dist)
{
   double v[3];

   novas_radec2vector (*ra, *dec, *dist, v);
   (void) push_doubles_as_array (v, 3);
}

static double interp_rad_vel (void)
{
   double pos[3], vel[3], vel_obs[3], d_obs[3], rv;
   Interp_Novas_object_t obj;

   if ((-1 == pop_doubles (d_obs, 3))
       || (-1 == pop_doubles_from_array (vel_obs, 3))
       || (-1 == pop_doubles_from_array (vel, 3))
       || (-1 == pop_doubles_from_array (pos, 3))
       || (-1 == pop_novas_object (&obj)))
     return -1.0;

   novas_rad_vel (&obj.object, pos, vel, vel_obs, d_obs[0], d_obs[1], d_obs[2], &rv);
   free_novas_object (&obj);
   return rv;
}

static double interp_refract (void)
{
   double zd, r;
   Novas_on_surface_t surface;
   short ref;

   if ((-1 == SLang_pop_double (&zd))
       || (-1 == SLang_pop_short (&ref))
       || (-1 == pop_novas_on_surface (&surface)))
     return -1.0;

   r = novas_refract (&surface, ref, zd);
   return r;
}

static void interp_spin (void)
{
   double ang, pos1[3], pos2[3];

   if ((-1 == pop_doubles_from_array (pos1, 3))
       || (-1 == SLang_pop_double (&ang)))
     return;

   novas_spin (ang, pos1, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

static void interp_starvectors (void)
{
   Interp_Novas_cat_entry_t cat;
   double pos[3], vel[3];

   if (-1 == pop_novas_cat_entry (&cat))
     return;

   novas_starvectors (&cat.cat_entry, pos, vel);
   (void) push_doubles_as_array (pos, 3);
   (void) push_doubles_as_array (vel, 3);
   free_novas_cat_entry (&cat);
}

static void interp_tdb2tt (double *tbd)
{
   double a[2];

   novas_tdb2tt (*tbd, a, a+1);
   (void) push_doubles (a, 2);
}

static void interp_terra (void)
{
   Novas_on_surface_t loc;
   double st, pos[3], vel[3];

   if ((-1 == SLang_pop_double (&st))
       || (-1 == pop_novas_on_surface (&loc)))
     return;

   novas_terra (&loc, st, pos, vel);
   (void) push_doubles_as_array (pos, 3);
   (void) push_doubles_as_array (vel, 3);
}

static void interp_transform_cat (void)
{
   double date_incat, date_newcat;
   Interp_Novas_cat_entry_t incat, newcat;
   char *new_cat_id;
   int ret;
   short opt;

   if (-1 == SLang_pop_slstring (&new_cat_id))
     return;
   if ((-1 == SLang_pop_double (&date_newcat))
       || (-1 == pop_novas_cat_entry (&incat)))
     {
	SLang_free_slstring (new_cat_id);
	return;
     }
   if ((-1 == SLang_pop_double (&date_incat))
       || (-1 == SLang_pop_short (&opt)))
     {
	free_novas_cat_entry (&incat);
	SLang_free_slstring (new_cat_id);
	return;
     }

   ret = novas_transform_cat (opt, date_incat, &incat.cat_entry,
			      date_newcat, new_cat_id, &newcat.cat_entry);
   free_novas_cat_entry (&incat);
   SLang_free_slstring (new_cat_id);

   if (ret == 0)
     {
	(void) push_novas_cat_entry (&newcat);
	return;
     }
   throw_error ("transform_cat", ret);
}

static void interp_transform_hip (void)
{
   Interp_Novas_cat_entry_t hip, hip2000;

   if (-1 == pop_novas_cat_entry (&hip))
     return;

   novas_transform_hip (&hip.cat_entry, &hip2000.cat_entry);
   free_novas_cat_entry (&hip);
   (void) push_novas_cat_entry (&hip2000);
}

static void interp_wobble (void)
{
   double tjd, xpyp[2], pos1[3], pos2[3];
   short dir;

   if ((-1 == pop_doubles_from_array (pos1, 3))
       || (-1 == pop_doubles (xpyp, 2))
       || (-1 == SLang_pop_short (&dir))
       || (-1 == SLang_pop_double (&tjd)))
     return;

   novas_wobble (tjd, dir, xpyp[0], xpyp[1], pos1, pos2);
   (void) push_doubles_as_array (pos2, 3);
}

#define V SLANG_VOID_TYPE
#define S SLANG_SHORT_TYPE
#define D SLANG_DOUBLE_TYPE
static SLang_Intrin_Fun_Type Module_Intrinsics [] =
{
   MAKE_INTRINSIC_0("app_star", interp_app_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("virtual_star", interp_virtual_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("astro_star", interp_astro_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("app_planet", interp_app_planet, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("virtual_planet", interp_virtual_planet, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("astro_planet", interp_astro_planet, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("topo_star", interp_topo_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("local_star", interp_local_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("topo_planet", interp_topo_planet, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("local_planet", interp_local_planet, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("mean_star", interp_mean_star, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("place", interp_place, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("julian_date", interp_julian_date, SLANG_DOUBLE_TYPE),
   MAKE_INTRINSIC_2("equ2gal", interp_equ2gal, V, D, D),
   MAKE_INTRINSIC_5("equ2ecl", interp_equ2ecl, V, D, S, S, D, D),
   MAKE_INTRINSIC_0("equ2hor", interp_equ2hor, V),
   MAKE_INTRINSIC_6("sidereal_time", interp_sidereal_time, V, D, D, D, S, S, S),
   MAKE_INTRINSIC_2("era", interp_era, V, D, D),
   MAKE_INTRINSIC_0("ephemeris", interp_ephemeris, V),
   MAKE_INTRINSIC_0("equ2ecl_vec", interp_equ2ecl_vec, V),
   MAKE_INTRINSIC_0("ecl2equ_vec", interp_ecl2equ_vec, V),
   MAKE_INTRINSIC_5("gcrs2equ", interp_gcrs2equ, V, D, S, S, D, D),
   MAKE_INTRINSIC_0("vector2radec", interp_vector2radec, V),
   MAKE_INTRINSIC_0("ter2cel", interp_ter2cel, V),
   MAKE_INTRINSIC_0("cel2ter", interp_cel2ter, V),
   MAKE_INTRINSIC_S("ephem_open", interp_ephem_open, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("ephem_close", interp_ephem_close, SLANG_VOID_TYPE),
   MAKE_INTRINSIC_0("aberration", interp_aberration, V),
   MAKE_INTRINSIC_0("bary2obs", interp_bary2obs, V),
   MAKE_INTRINSIC_1("cal_date", interp_cal_date, V, D),
   MAKE_INTRINSIC_4("cel_pole", interp_cel_pole, V, D, S, D, D),
   MAKE_INTRINSIC_2("cio_array", interp_cio_array, V, D, SLANG_LONG_TYPE),
   MAKE_INTRINSIC_4("cio_basis", interp_cio_basis, V, D, D, S, S),
   MAKE_INTRINSIC_2("cio_location", interp_cio_location, V, D, S),
   MAKE_INTRINSIC_2("cio_ra", interp_cio_ra, D, D, S),
   MAKE_INTRINSIC_0("d_light", interp_d_light, D),
   MAKE_INTRINSIC_3("ee_ct", interp_ee_ct, D, D, D, S),
   MAKE_INTRINSIC_2("e_tilt", interp_e_tilt, V, D, S),
   MAKE_INTRINSIC_0("frame_tie", interp_frame_tie, V),
   MAKE_INTRINSIC_1("fund_args", interp_fund_args, V, D),
   MAKE_INTRINSIC_0("geo_posvel", interp_geo_posvel, V),
   MAKE_INTRINSIC_0("grav_def", interp_grav_def, V),
   MAKE_INTRINSIC_0("grav_vec", interp_grav_vec, V),
   MAKE_INTRINSIC_3("ira_equinox", interp_ira_equinox, D, D, S, S),
   MAKE_INTRINSIC_0("light_time", interp_light_time, V),
   MAKE_INTRINSIC_0("limb_angle", interp_limb_angle, V),
   MAKE_INTRINSIC_1("mean_obliq", interp_mean_obliq, D, D),
   MAKE_INTRINSIC_1("norm_ang", interp_norm_ang, D, D),
   MAKE_INTRINSIC_0("nutation", interp_nutation, V),
   MAKE_INTRINSIC_2("nutation_angles", interp_nutation_angles, V, D, S),
   MAKE_INTRINSIC_0("planet_ephemeris", interp_planet_ephemeris, V),
   MAKE_INTRINSIC_0("precession", interp_precession, V),
   MAKE_INTRINSIC_0("proper_motion", interp_proper_motion, V),
   MAKE_INTRINSIC_3("radec2vector", interp_radec2vector, V, D, D, D),
   MAKE_INTRINSIC_0("rad_vel", interp_rad_vel, D),
   MAKE_INTRINSIC_0("refract", interp_refract, D),
   MAKE_INTRINSIC_0("spin", interp_spin, V),
   MAKE_INTRINSIC_0("starvectors", interp_starvectors, V),
   MAKE_INTRINSIC_1("tdb2tt", interp_tdb2tt, V, D),
   MAKE_INTRINSIC_0("terra", interp_terra, V),
   MAKE_INTRINSIC_0("transform_cat", interp_transform_cat, V),
   MAKE_INTRINSIC_0("transform_hip", interp_transform_hip, V),
   MAKE_INTRINSIC_0("wobble", interp_wobble, V),

   MAKE_INTRINSIC_0("novas_set_throw_callback", interp_set_throw_callback, V),
   SLANG_END_INTRIN_FUN_TABLE
};

#undef V
#undef S
#undef D
static SLang_Intrin_Var_Type Module_Variables [] =
{
   MAKE_VARIABLE("_slnovas_module_version_string", &Module_Version_String, SLANG_STRING_TYPE, 1),
   SLANG_END_INTRIN_VAR_TABLE
};

static SLang_IConstant_Type Module_IConstants [] =
{
   MAKE_ICONSTANT("_slnovas_module_version", MODULE_VERSION_NUMBER),
   SLANG_END_ICONST_TABLE
};

static SLang_DConstant_Type Module_DConstants [] =
{
   SLANG_END_DCONST_TABLE
};

int init_slnovas_module_ns (char *ns_name)
{
   SLang_NameSpace_Type *ns = SLns_create_namespace (ns_name);
   if (ns == NULL)
     return -1;

   if (
       (-1 == SLns_add_intrin_var_table (ns, Module_Variables, NULL))
       || (-1 == SLns_add_intrin_fun_table (ns, Module_Intrinsics, NULL))
       || (-1 == SLns_add_iconstant_table (ns, Module_IConstants, NULL))
       || (-1 == SLns_add_dconstant_table (ns, Module_DConstants, NULL))
       )
     return -1;

   return 0;
}

/* This function is optional */
void deinit_slnovas_module (void)
{
}
