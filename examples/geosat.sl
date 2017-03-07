require ("slnovas");

private variable GEO_SAT_RADIUS = 42163.968; % km
private variable KM_PER_AU = 149597871.0;

private define lon_lat_to_vec (lon, lat)
{
   variable torad = PI/180.0;
   lon *= torad; lat *= torad;
   variable cos_lat = cos(lat);
   return [cos_lat*cos(lon), cos_lat*sin(lon), sin(lat)];
}

private define angle_between_vectors (a, b)
{
   a /= hypot(a);
   b /= hypot(b);
   variable cos_theta = sum(a*b);
   if (-0.7 <= cos_theta <= 0.7)
     return acos (cos_theta);
   variable c =
     hypot ([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]);

   return atan2 (c, cos_theta);
}

define slsh_main ()
{
   variable ephemfile = "../libnovas/JPLEPH";
   if (NULL == stat_file (ephemfile))
     ephemfile = "JPLEPH";

   (,,) = ephem_open (ephemfile);

   variable boresight_lon = -100.0;    %  geodetic
   variable boresight_lat = 36.5;

   variable sat_lon = boresight_lon;   %  geodetic
   variable sat_lat = 0.0;

   % Here, sat_pos is the position of a geo-synchronous satellite that
   % is stationary with respect to the ITRS system.
   variable sat_pos = GEO_SAT_RADIUS * lon_lat_to_vec(sat_lon, sat_lat);
   variable bs_obs = make_observer_on_surface (boresight_lat, boresight_lon, 0, 0, 0);
   variable geocenter = make_observer_at_geocenter ();
   variable sun = make_object (0, 10, "sun", NULL);

   variable accuracy = 0;	       %  full

   variable year0 = 2016, month0 = 3, day0 = 21;
   variable year1 = 2016, month1 = 3, day1 = 22;

   variable jd_utc0 = julian_date (year0, month0, day0, 0.0); % days
   variable jd_utc1 = julian_date (year1, month1, day1, 0.0); % days

   variable jd_utc = jd_utc0;

   while (jd_utc < jd_utc1)
     {
	variable leap_secs = 36;	       %  current
	variable ut1_utc = -0.21325;      %  predicted UT1-UTC

	variable jd_tt = jd_utc + (leap_secs + 32.184)/86400;
	variable jd_ut1 = jd_utc + ut1_utc/86400;
	variable delta_t = 32.184 + leap_secs - ut1_utc;

	% Convert to GCRS
	variable sat_gcrs, bs_gcrs, sun_gcrs, sun_place, angle;
	sat_gcrs = ter2cel (jd_ut1, 0.0, delta_t, 1, accuracy, 0, 0.0, 0.0, sat_pos);
	(bs_gcrs, ) = geo_posvel (jd_tt, delta_t, accuracy, bs_obs);
	sun_place = place (jd_tt, sun, geocenter, delta_t, 0, accuracy);
	sun_gcrs = sun_place.dis*sun_place.r_hat * KM_PER_AU;

	angle = angle_between_vectors (bs_gcrs-sat_gcrs, sun_gcrs-sat_gcrs);

	variable cal = cal_date (jd_utc);
	variable hour = int(cal.hour);
	variable min = int (60*(cal.hour - hour));
	variable sec = int (3600*(cal.hour - hour - min/60.0));
	() = fprintf (stdout, "%S %S %d-%02d-%02dT%02d:%02d:%02d\n",
		      cal.hour,
		      180/PI*angle,
		      cal.year, cal.month, cal.day,
		      hour, min, sec);


	jd_utc += 1*60/86400.0;
     }
}
