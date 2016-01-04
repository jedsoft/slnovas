prepend_to_slang_load_path ("..");
set_import_module_path ("..:" + get_import_module_path ());

require ("slnovas.sl");
private variable N_STARS = 3, N_TIMES = 4;

define slsh_main ()
{
   %Main function to check out many parts of NOVAS-C by calling
   %function 'topo_star' with version 1 of function 'solarsystem'.

   variable error = 0;
   variable accuracy = 0;
   variable i, j, de_num;

   % 'deltat' is the difference in time scales, TT - UT1.
   % The array 'tjd' contains four selected Julian dates at which the
   % star positions will be evaluated.


   variable deltat = 60.0;
   variable tjd = [2450203.5, 2450203.5, 2450417.5, 2450300.5];
   variable jd_beg, jd_end, ra, dec;

   % Hipparcos (ICRS) catalog data for three selected stars.

   variable stars =
     {
	make_cat_entry ("POLARIS", "HIP",   0,  2.530301028,  89.264109444,
			44.22, -11.75,  7.56, -17.4),
	make_cat_entry ("Delta ORI", "HIP", 1,  5.533444639,  -0.299091944,
			1.67,   0.56,  3.56,  16.0),
	make_cat_entry ("Theta CAR", "HIP", 2, 10.715944806, -64.394450000,
			-18.87, 12.06,  7.43,  24.0),
     };

   % The observer's terrestrial coordinates (latitude, longitude, height).
   variable geo_loc = make_on_surface (45.0, -75.0, 0.0, 10.0, 1010.0);

   % Open the JPL ephemeris file.

   % Open it twice to check for leaks
   (jd_beg, jd_end, de_num) = ephem_open ("../../libnovas/JPLEPH");
   ephem_close ();
   (jd_beg, jd_end, de_num) = ephem_open ("../../libnovas/JPLEPH");

   ()=printf ("JPL Ephemeris DE%d open. jd_beg = %10.2f  jd_end = %10.2f\n",
	      de_num, jd_beg, jd_end);
   ()=printf ("\n");

   % Compute the topocentric places of the three stars at the four
   % selected Julian dates.

   for (i = 0; i < N_TIMES; i++)
   {
      for (j = 0; j < N_STARS; j++)
      {
	 (ra, dec) = topo_star (tjd[i],deltat,stars[j],geo_loc, accuracy);
	 () = printf ("JD = %f  Star = %s\n", tjd[i], stars[j].starname);
	 () = printf ("RA = %12.9f  Dec = %12.8f\n", ra, dec);
	 () = printf ("\n");
      }
      () = printf ("\n");
   }

   ephem_close ();
}
