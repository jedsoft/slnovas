/* This file was derived from Cdist/novas.h */

#ifndef LIBNOVAS_H
# define LIBNOVAS_H

/*
 struct Novas_cat_entry_t:  basic astrometric data for any celestial object
 located outside the solar system; the catalog
 data for a star
 *
 starname[NOVAS_SIZE_OF_OBJ_NAME] = name of celestial object
 catalog[NOVAS_SIZE_OF_CAT_NAME]  = catalog designator (e.g., HIP)
 starnumber                 = integer identifier assigned to object
 ra                         = ICRS right ascension (hours)
 dec                        = ICRS declination (degrees)
 promora                    = ICRS proper motion in right ascension
 (milliarcseconds/year)
 promodec                   = ICRS proper motion in declination
 (milliarcseconds/year)
 parallax                   = parallax (milliarcseconds)
 radialvelocity             = radial velocity (km/s)
 *
 NOVAS_SIZE_OF_OBJ_NAME and NOVAS_SIZE_OF_CAT_NAME are defined below.  Each is the
 number of characters in the string (the string length) plus the null
 terminator.
 */

# define NOVAS_SIZE_OF_OBJ_NAME 51
# define NOVAS_SIZE_OF_CAT_NAME 4

typedef struct
{
   char starname[NOVAS_SIZE_OF_OBJ_NAME];
   char catalog[NOVAS_SIZE_OF_CAT_NAME];
   long int starnumber;
   double ra;
   double dec;
   double promora;
   double promodec;
   double parallax;
   double radialvelocity;
} Novas_cat_entry_t;

/*
 Novas_object_t:    specifies the celestial object of interest
 *
 type              = type of object
 = 0 ... major planet, Pluto, Sun, or Moon
 = 1 ... minor planet
 = 2 ... object located outside the solar system
 (star, nebula, galaxy, etc.)
 number            = object number
 For 'type' = 0: Mercury = 1, ..., Pluto = 9,
 Sun = 10, Moon = 11
 For 'type' = 1: minor planet number
 For 'type' = 2: set to 0 (object is
 fully specified in 'struct cat_entry')
 name              = name of the object (limited to
 (NOVAS_SIZE_OF_OBJ_NAME - 1) characters)
 star              = basic astrometric data for any celestial object
 located outside the solar system; the catalog
 data for a star
 */

typedef struct
{
   short int type;
   short int number;
   char name[NOVAS_SIZE_OF_OBJ_NAME];
   Novas_cat_entry_t star;
} Novas_object_t;

/*
 Novas_on_surface_t: data for an observer's location on the surface of
 the Earth.  The atmospheric parameters are used
 only by the refraction function called from
 function 'equ2hor'. Additional parameters can be
 added to this structure if a more sophisticated
 refraction model is employed.
 *
 latitude           = geodetic (ITRS) latitude; north positive (degrees)
 longitude          = geodetic (ITRS) longitude; east positive (degrees)
 height             = height of the observer (meters)
 temperature        = temperature (degrees Celsius)
 pressure           = atmospheric pressure (millibars)
 */

typedef struct
{
   double latitude;
   double longitude;
   double height;
   double temperature;
   double pressure;
} Novas_on_surface_t;

/*
 Novas_in_space_t:   data for an observer's location on a near-Earth
 spacecraft
 *
 sc_pos[3]          = geocentric position vector (x, y, z), components
 in km
 sc_vel[3]          = geocentric velocity vector (x_dot, y_dot,
 z_dot), components in km/s
 *
 Both vectors with respect to true equator and
 equinox of date
 */

typedef struct
{
   double sc_pos[3];
   double sc_vel[3];
} Novas_in_space_t;

/*
 Novas_observer_t:   data specifying the location of the observer
 *
 where              = integer code specifying location of observer
 = 0: observer at geocenter
 = 1: observer on surface of earth
 = 2: observer on near-earth spacecraft
 on_surface         = structure containing data for an observer's
 location on the surface of the Earth (where = 1)
 near_earth         = data for an observer's location on a near-Earth
 spacecraft (where = 2)
 */

typedef struct
{
   short int where;
   Novas_on_surface_t on_surf;
   Novas_in_space_t near_earth;
} Novas_observer_t;

/*
 Novas_sky_pos_t:    data specifying a celestial object's place on the
 sky; contains the output from function 'place'
 *
 r_hat[3]           = unit vector toward object (dimensionless)
 ra                 = apparent, topocentric, or astrometric
 right ascension (hours)
 dec                = apparent, topocentric, or astrometric
 declination (degrees)
 dis                = true (geometric, Euclidian) distance to solar
 system body or 0.0 for star (AU)
 rv                 = radial velocity (km/s)
 */

typedef struct
{
   double r_hat[3];
   double ra;
   double dec;
   double dis;
   double rv;
} Novas_sky_pos_t;

/*
 Novas_ra_of_cio_t:  right ascension of the Celestial Intermediate
 Origin (CIO) with respect to the GCRS
 *
 jd_tdb             = TDB Julian date
 ra_cio             = right ascension of the CIO with respect
 to the GCRS (arcseconds)
 */

typedef struct
{
   double jd_tdb;
   double ra_cio;
} Novas_ra_of_cio_t;

/*
 Define "origin" constants.
 */

# define NOVAS_BARYC  0
# define NOVAS_HELIOC 1

extern short int novas_app_star (double jd_tt, Novas_cat_entry_t *star,
				 short int accuracy,
				 double *ra, double *dec);

extern short int novas_virtual_star (double jd_tt, Novas_cat_entry_t *star,
				     short int accuracy,
				     double *ra, double *dec);

extern short int novas_astro_star (double jd_tt, Novas_cat_entry_t *star,
				   short int accuracy,
				   double *ra, double *dec);

extern short int novas_app_planet (double jd_tt, Novas_object_t *ss_body,
				   short int accuracy,
				   double *ra, double *dec, double *dis);

extern short int novas_virtual_planet (double jd_tt, Novas_object_t *ss_body,
				       short int accuracy,
				       double *ra, double *dec, double *dis);

extern short int novas_astro_planet (double jd_tt, Novas_object_t *ss_body,
				     short int accuracy,
				     double *ra, double *dec, double *dis);

extern short int novas_topo_star (double jd_tt, double delta_t,
				  Novas_cat_entry_t *star,
				  Novas_on_surface_t *position, short int accuracy,
				  double *ra, double *dec);

extern short int novas_local_star (double jd_tt, double delta_t,
				   Novas_cat_entry_t *star,
				   Novas_on_surface_t *position, short int accuracy,
				   double *ra, double *dec);

extern short int novas_topo_planet (double jd_tt, Novas_object_t *ss_body,
				    double delta_t,
				    Novas_on_surface_t *position, short int accuracy,
				    double *ra, double *dec, double *dis);

extern short int novas_local_planet (double jd_tt, Novas_object_t *ss_body,
				     double delta_t, Novas_on_surface_t *position,
				     short int accuracy,
				     double *ra, double *dec, double *dis);

extern short int novas_mean_star (double jd_tt, double ra, double dec,
				  short int accuracy,
				  double *ira, double *idec);

extern short int novas_place (double jd_tt, Novas_object_t *cel_object,
			      Novas_observer_t *location, double delta_t,
			      short int coord_sys, short int accuracy,
			      Novas_sky_pos_t *output);

extern void novas_equ2gal (double rai, double deci,
			   double *glon, double *glat);

extern short int novas_equ2ecl (double jd_tt, short int coord_sys,
				short int accuracy, double ra, double dec,
				double *elon, double *elat);

extern short int novas_equ2ecl_vec (double jd_tt, short int coord_sys,
				    short int accuracy, double *pos1,
				    double *pos2);

extern short int novas_ecl2equ_vec (double jd_tt, short int coord_sys,
				    short int accuracy, double *pos1,
				    double *pos2);

extern void novas_equ2hor (double jd_ut1, double delta_t, short int accuracy,
			   double xp, double yp, Novas_on_surface_t *location, double ra,
			   double dec, short int ref_option,
			   double *zd, double *az, double *rar, double *decr);

extern short int novas_gcrs2equ (double jd_tt, short int coord_sys,
				 short int accuracy, double rag, double decg,
				 double *ra, double *dec);

extern short int novas_sidereal_time (double jd_high, double jd_low,
				      double delta_t, short int gst_type,
				      short int method, short int accuracy,
				      double *gst);

extern double novas_era (double jd_high, double jd_low);

extern short int novas_ter2cel (double jd_ut_high, double jd_ut_low,
				double delta_t, short int method,
				short int accuracy, short int option, double xp,
				double yp, double *vec1,
				double *vec2);

extern short int novas_cel2ter (double jd_ut_high, double jd_ut_low,
				double delta_t, short int method,
				short int accuracy, short int option,
				double xp, double yp, double *vec1,
				double *vec2);

extern void novas_spin (double angle, double *pos1,
			double *pos2);

extern void novas_wobble (double tjd, short int direction, double xp, double yp,
			  double *pos1,
			  double *pos2);

extern void novas_terra (Novas_on_surface_t *location, double st,
			 double *pos, double *vel);

extern void novas_e_tilt (double jd_tdb, short int accuracy,
			  double *mobl, double *tobl, double *ee, double *dpsi,
			  double *deps);

extern short int novas_cel_pole (double tjd, short int type, double dpole1,
				 double dpole2);

extern double novas_ee_ct (double jd_high, double jd_low, short int accuracy);

extern void novas_frame_tie (double *pos1, short int direction,

			     double *pos2);

extern void novas_proper_motion (double jd_tdb1, double *pos, double *vel,
				 double jd_tdb2,

				 double *pos2);

extern void novas_bary2obs (double *pos, double *pos_obs,

			    double *pos2, double *lighttime);

extern short int novas_geo_posvel (double jd_tt, double delta_t,
				   short int accuracy, Novas_observer_t *obs,
				   double *pos, double *vel);

extern short int novas_light_time (double jd_tdb, Novas_object_t *ss_object,
				   double pos_obs[3], double tlight0,
				   short int accuracy,
				   double pos[3], double *tlight);

extern double novas_d_light (double *pos1, double *pos_obs);

extern short int novas_grav_def (double jd_tdb, short int loc_code,
				 short int accuracy, double *pos1, double *pos_obs,
				 double *pos2);

extern void novas_grav_vec (double *pos1, double *pos_obs, double *pos_body,
			    double rmass,
			    double *pos2);

extern void novas_aberration (double *pos, double *ve, double lighttime,
			      double *pos2);

extern void novas_rad_vel (Novas_object_t *cel_object, double *pos, double *vel,
			   double *vel_obs, double d_obs_geo, double d_obs_sun,
			   double d_obj_sun,
			   double *rv);

extern short int novas_precession (double jd_tdb1, double *pos1, double jd_tdb2,
				   double *pos2);

extern void novas_nutation (double jd_tdb, short int direction, short int accuracy,
			    double *pos,
			    double *pos2);

extern void novas_nutation_angles (double t, short int accuracy,
				   double *dpsi, double *deps);

extern void novas_fund_args (double t,
			     double a[5]);

extern double novas_mean_obliq (double jd_tdb);

extern short int novas_vector2radec (double *pos,
				     double *ra, double *dec);

extern void novas_radec2vector (double ra, double dec, double dist,
				double *vector);

extern void novas_starvectors (Novas_cat_entry_t *star,
			       double *pos, double *vel);

extern void novas_tdb2tt (double tdb_jd,
			  double *tt_jd, double *secdiff);

extern short int novas_cio_ra (double jd_tt, short int accuracy,
			       double *ra_cio);

extern short int novas_cio_location (double jd_tdb, short int accuracy,
				     double *ra_cio, short int *ref_sys);

extern short int novas_cio_basis (double jd_tdb, double ra_cio, short int ref_sys,
				  short int accuracy,
				  double *x, double *y, double *z);

extern short int novas_cio_array (double jd_tdb, long int n_pts,
				  Novas_ra_of_cio_t *cio);

extern double novas_ira_equinox (double jd_tdb, short int equinox,
				 short int accuracy);

extern short int novas_ephemeris (double jd[2], Novas_object_t *cel_obj, short int origin,
				  short int accuracy,
				  double *pos, double *vel);

extern void novas_transform_hip (Novas_cat_entry_t *hipparcos,
				 Novas_cat_entry_t *hip_2000);

extern short int novas_transform_cat (short int option, double date_incat,
				      Novas_cat_entry_t *incat, double date_newcat,
				      char newcat_id[NOVAS_SIZE_OF_CAT_NAME],
				      Novas_cat_entry_t *newcat);

extern void novas_limb_angle (double pos_obj[3], double pos_obs[3],
			      double *limb_ang, double *nadir_ang);

extern double novas_refract (Novas_on_surface_t *location, short int ref_option,
			     double zd_obs);

extern double novas_julian_date (short int year, short int month, short int day,
				 double hour);

extern void novas_cal_date (double tjd,
			    short int *year, short int *month, short int *day,
			    double *hour);

extern double novas_norm_ang (double angle);

extern short int novas_make_cat_entry (char star_name[NOVAS_SIZE_OF_OBJ_NAME],
				       char catalog[NOVAS_SIZE_OF_CAT_NAME],
				       long int star_num, double ra, double dec,
				       double pm_ra, double pm_dec, double parallax,
				       double rad_vel,
				       Novas_cat_entry_t *star);

extern short int novas_make_object (short int type, short int number,
				    char name[NOVAS_SIZE_OF_OBJ_NAME],
				    Novas_cat_entry_t *star_data,
				    Novas_object_t *cel_obj);

extern short int novas_make_observer (short int where, Novas_on_surface_t *obs_surface,
				      Novas_in_space_t *obs_space,
				      Novas_observer_t *obs);

extern void novas_make_observer_at_geocenter (
					      Novas_observer_t *obs_at_geocenter);

extern void novas_make_observer_on_surface (double latitude, double longitude,
					    double height, double temperature,
					    double pressure,
					    Novas_observer_t *obs_on_surface);

extern void novas_make_observer_in_space (double sc_pos[3], double sc_vel[3],
					  Novas_observer_t *obs_in_space);

extern void novas_make_on_surface (double latitude, double longitude,
				   double height,
				   double temperature, double pressure,
				   Novas_on_surface_t *obs_surface);

extern void novas_make_in_space (double sc_pos[3], double sc_vel[3],
				 Novas_in_space_t *obs_space);

extern short int novas_ephem_open (char *ephem_name,
			    double *jd_begin, double *jd_end,
			    short int *de_number);
extern short int novas_ephem_close (void);

extern short int novas_planet_ephemeris (double tjd[2], short int target,
				  short int center,
				  double *position, double *velocity);

#endif				       /* LIBNOVAS_H */
