#ifndef P_LIBNOVAS_H
#define P_LIBNOVAS_H

/* novas.h */
#define _NOVAS_	1		       /* prevent inclusion of novas.h */
#define SIZE_OF_OBJ_NAME NOVAS_SIZE_OF_OBJ_NAME
#define SIZE_OF_CAT_NAME NOVAS_SIZE_OF_CAT_NAME

#define cat_entry Novas_cat_entry_t
#define object Novas_object_t
#define on_surface Novas_on_surface_t
#define in_space Novas_in_space_t
#define observer Novas_observer_t
#define sky_pos Novas_sky_pos_t
#define ra_of_cio Novas_ra_of_cio_t

#define BARYC NOVAS_BARYC
#define HELIOC NOVAS_HELIOC

#ifdef LIBNOVAS_SOURCE
static double *readeph (int mp, char *name, double jd, int *err);
#endif

#define app_star novas_app_star
#define virtual_star novas_virtual_star
#define astro_star novas_astro_star
#define app_planet novas_app_planet
#define virtual_planet novas_virtual_planet
#define astro_planet novas_astro_planet
#define topo_star novas_topo_star
#define local_star novas_local_star
#define topo_planet novas_topo_planet
#define local_planet novas_local_planet
#define mean_star novas_mean_star
#define place novas_place
#define equ2gal novas_equ2gal
#define equ2ecl novas_equ2ecl
#define equ2ecl_vec novas_equ2ecl_vec
#define ecl2equ_vec novas_ecl2equ_vec
#define equ2hor novas_equ2hor
#define gcrs2equ novas_gcrs2equ
#define sidereal_time novas_sidereal_time
#define era novas_era
#define ter2cel novas_ter2cel
#define cel2ter novas_cel2ter
#define spin novas_spin
#define wobble novas_wobble
#define terra novas_terra
#define e_tilt novas_e_tilt
#define cel_pole novas_cel_pole
#define ee_ct novas_ee_ct
#define frame_tie novas_frame_tie
#define proper_motion novas_proper_motion
#define bary2obs novas_bary2obs
#define geo_posvel novas_geo_posvel
#define light_time novas_light_time
#define d_light novas_d_light
#define grav_def novas_grav_def
#define grav_vec novas_grav_vec
#define aberration novas_aberration
#define rad_vel novas_rad_vel
#define precession novas_precession
#define nutation novas_nutation
#define nutation_angles novas_nutation_angles
#define fund_args novas_fund_args
#define vector2radec novas_vector2radec
#define radec2vector novas_radec2vector
#define starvectors novas_starvectors
#define tdb2tt novas_tdb2tt
#define cio_ra novas_cio_ra
#define cio_location novas_cio_location
#define cio_basis novas_cio_basis
#define cio_array novas_cio_array
#define ira_equinox novas_ira_equinox
#define ephemeris novas_ephemeris
#define transform_hip novas_transform_hip
#define transform_cat novas_transform_cat
#define limb_angle novas_limb_angle
#define refract novas_refract
#define julian_date novas_julian_date
#define cal_date novas_cal_date
#define norm_ang novas_norm_ang
#define make_cat_entry novas_make_cat_entry
#define make_object novas_make_object
#define make_observer novas_make_observer
#define make_observer_at_geocenter novas_make_observer_at_geocenter
#define make_observer_on_surface novas_make_observer_on_surface
#define make_observer_in_space novas_make_observer_in_space
#define make_on_surface novas_make_on_surface
#define make_in_space novas_make_in_space

/* solarsystem.h */
#define _SOLSYS_ 1		       /* prevent inclusion of solarsystem.h */
#ifdef LIBNOVAS_SOURCE
static short int solarsystem (double tjd, short int body, short int origin,
			      double *position, double *velocity);

static short int solarsystem_hp (double tjd[2], short body, short origin,
				 double *position, double *velocity);
#endif

/* novacon.h */
#define _CONSTS_ 1		       /* prevent inclusion */
#define T0 NOVAS_CONST_T0
#define FN1 NOVAS_CONST_FN1
#define FN0 NOVAS_CONST_FN0
#define T0 NOVAS_CONST_T0
#define C NOVAS_CONST_C
#define AU_SEC NOVAS_CONST_AU_SEC
#define C_AUDAY NOVAS_CONST_C_AUDAY
#define AU NOVAS_CONST_AU
#define AU_KM NOVAS_CONST_AU_KM
#define GS NOVAS_CONST_GS
#define GE NOVAS_CONST_GE
#define ERAD NOVAS_CONST_ERAD
#define F NOVAS_CONST_F
#define ANGVEL NOVAS_CONST_ANGVEL
#define RMASS NOVAS_CONST_RMASS
#define TWOPI NOVAS_CONST_TWOPI
#define ASEC360 NOVAS_CONST_ASEC360
#define ASEC2RAD NOVAS_CONST_ASEC2RAD
#define DEG2RAD NOVAS_CONST_DEG2RAD
#define RAD2DEG NOVAS_CONST_RAD2DEG

/* nutation.h */
#define _NUTATION_ 1
#ifdef LIBNOVAS_SOURCE
static void iau2000a (double jd_high, double jd_low,
		      double *dpsi, double *deps);

static void iau2000b (double jd_high, double jd_low,
		      double *dpsi, double *deps);

static void nu2000k (double jd_high, double jd_low,
		      double *dpsi, double *deps);
#endif

/* eph_manager.h */
#if (NOVAS_SOLSYS_MODEL == 1)
# define _EPHMAN_  1

# define KM _pNOVAS_EPHMAN_KM
# define IPT _pNOVAS_EPHMAN_IPT
# define LPT _pNOVAS_EPHMAN_LPT
# define NRL _pNOVAS_EPHMAN_NRL
# define NP _pNOVAS_EPHMAN_NP
# define NV _pNOVAS_EPHMAN_NV
# define RECORD_LENGTH _pNOVAS_EPHMAN_RECORD_LENGTH
# define SS _pNOVAS_EPHMAN_SS
# define JPLAU _pNOVAS_EPHMAN_JPLAU
# define PC _pNOVAS_EPHMAN_PC
# define VC _pNOVAS_EPHMAN_VC
# define TWOT _pNOVAS_EPHMAN_TWOT
# define EM_RATIO _pNOVAS_EPHMAN_EM_RATIO
# define BUFFER _pNOVAS_EPHMAN_BUFFER
# define EPHFILE _pNOVAS_EPHMAN_EPHFILE

#define planet_ephemeris novas_planet_ephemeris

#ifdef LIBNOVAS_SOURCE
static short int ephem_open (char *ephem_name,
			    double *jd_begin, double *jd_end,
			    short int *de_number);
static short int ephem_close (void);
static short int state (double *jed, short int target,
                 double *target_pos, double *target_vel);
static void interpolate (double *buf, double *t, long int ncm, long int na,
                  double *position, double *velocity);
static void split (double tt, double *fr);
#endif

#endif				       /* NOVAS_SOLSYS_MODEL==1 */

#endif 				       /* P_LIBNOVAS_H */
