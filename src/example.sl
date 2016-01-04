require ("slnovas");

private variable T0 = 2451545.00000000;
private variable DEG2RAD = 0.017453292519943296;

public define slsh_main ()
{
   variable year = 2015, month = 12, day = 15, leap_secs = 33, accuracy = 0;

   variable hour = 10 + 5.0/60 + 5;
   variable ut1_utc = -0.387845;

   variable latitude = 42.3605;
   variable longitude = -71.0596;
   variable height = 0.0;
   variable temperature = 10.0;
   variable pressure = 1010.0;

   variable x_pole = -0.002;
   variable y_pole = +0.529;

   variable geo_loc, obs_loc, star, dummy_star, mars, moon;

   geo_loc = make_on_surface (latitude, longitude, height, temperature, pressure);
   obs_loc = make_observer_on_surface (latitude, longitude, height, temperature, pressure);
   star = make_cat_entry ("GMB 1830","FK6",1307,11.88299133,37.71867646,
		   4003.27,-5815.07,109.21,-98.8);

   dummy_star = make_cat_entry ("DUMMY","xxx",0,0.0,0.0,0.0,0.0,0.0,0.0);
   moon = make_object (0, 11, "Moon", dummy_star);
   mars = make_object (0, 4, "Mars", dummy_star);

   variable jd_beg, jd_en, de_num;
   (jd_beg, jd_en, de_num) = ephem_open ("../libnovas/JPLEPH");

   ()=printf ("NOVAS Sample Calculations\n");
   ()=printf ("-------------------------\n");
   ()=printf ("\n");

   ()=printf ("Geodetic location:\n");
   ()=printf ("%15.10f        %15.10f        %15.10f\n\n", geo_loc.longitude,
	   geo_loc.latitude, geo_loc.height);

   variable jd_utc = julian_date (year,month,day,hour);

   variable jd_tt = jd_utc + (leap_secs + 32.184) / 86400.0;
   variable jd_ut1 = jd_utc + ut1_utc / 86400.0;
   variable delta_t = 32.184 + leap_secs - ut1_utc;

   variable jd_tdb = jd_tt;         % Approximation good to 0.0017 seconds.

   ()=printf ("TT and UT1 Julian Dates and Delta-T:\n");
   ()=printf ("%15.6f        %15.6f        %16.11f\n", jd_tt, jd_ut1, delta_t);
   ()=printf ("\n");

   variable ra, dec, rat, dect, dis, dist;

   (ra, dec) = app_star (jd_tt,star,accuracy);
   (rat, dect) = topo_star (jd_tt,delta_t,star,geo_loc, accuracy);

   ()=printf ("FK6 1307 geocentric and topocentric positions:\n");
   ()=printf ("%15.10f        %15.10f\n", ra, dec);
   ()=printf ("%15.10f        %15.10f\n", rat, dect);
   ()=printf ("\n");

   % Apparent and topocentric place of the Moon.

   (ra, dec, dis) = app_planet (jd_tt, moon, accuracy);

   (rat, dect, dist) = topo_planet (jd_tt,moon,delta_t,geo_loc,accuracy);

   ()=printf ("Moon geocentric and topocentric positions:\n");
   ()=printf ("%15.10f        %15.10f        %15.12f\n", ra, dec, dis);
   ()=printf ("%15.10f        %15.10f        %15.12f\n", rat, dect, dist);

   %
   % Topocentric place of the Moon using function 'place'-- should be
   % same as above.
   %

   variable t_place = place (jd_tt,moon,obs_loc,delta_t,1,accuracy);

   ()=printf ("%15.10f        %15.10f        %15.12f\n", t_place.ra,t_place.dec,
      t_place.dis);
   ()=printf ("\n");

   % Position of the Moon in local horizon coordinates.  (Polar motion
   % ignored here.)

   variable zd, az, rar, decr;
   (zd, az, rar, decr) = equ2hor (jd_ut1,delta_t,accuracy,0.0,0.0,
				  geo_loc,rat,dect,1);

   ()=printf ("Moon zenith distance and azimuth:\n");
   ()=printf ("%15.10f        %15.10f (decr=%.15f)\n", zd, az, decr);
   ()=printf ("\n");

   % Greenwich and local apparent sidereal time and Earth Rotation Angle.

   variable gast = sidereal_time (jd_ut1,0.0,delta_t,1,1,accuracy);

   variable last = gast + geo_loc.longitude / 15.0;
   if (last >= 24.0)
      last -= 24.0;
   if (last < 0.0)
      last += 24.0;

   variable theta = era (jd_ut1,0.0);

   ()=printf ("Greenwich and local sidereal time and Earth Rotation Angle:\n");
   ()=printf ("%16.11f        %16.11f        %15.10f\n", gast, last, theta);
   ()=printf ("\n");

   % Heliocentric position of Mars in BCRS.
   % TDB ~ TT approximation could lead to error of ~50 m in position of Mars.

   variable jd = [jd_tdb, 0.0];
   variable pos, vel;
   (pos, vel) = ephemeris (jd,mars,1,accuracy);

   variable pose;
   pose = equ2ecl_vec (T0,2,accuracy, pos);

   variable elat, elon;
   (elon, elat) = vector2radec (pose);
   elon *= 15.0;

   variable r = hypot(pose);

   ()=printf ("Mars heliocentric ecliptic longitude and latitude and radius vector:\n");
   ()=printf ("%15.10f        %15.10f        %15.12f\n", elon, elat, r);
   ()=printf ("\n");

   % Terrestrial to celestial transformation.

   variable lon_rad = geo_loc.longitude * DEG2RAD;
   variable lat_rad = geo_loc.latitude * DEG2RAD;
   variable sin_lon = sin (lon_rad);
   variable cos_lon = cos (lon_rad);
   variable sin_lat = sin (lat_rad);
   variable cos_lat = cos (lat_rad);

   % Form vector toward local zenith (orthogonal to ellipsoid) in ITRS.
   variable vter = [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat];

   % Transform vector to GCRS.

   variable vcel = ter2cel (jd_ut1,0.0,delta_t,1,accuracy,0,x_pole,y_pole,vter);

   (ra, dec) = vector2radec (vcel);

   ()=printf ("Direction of zenith vector (RA & Dec) in GCRS:\n");
   ()=printf ("%15.10f        %15.10f\n", ra, dec);
   ()=printf ("\n");

   ephem_close();
}

