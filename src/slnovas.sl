import ("slnovas");

% Error codes from the novas library are context-sensitive. The
% following code is designed to extract a meaningful error message from
% the called function context.
private variable Error_Info_Type = struct
{
   called,			       %  list of called functions and error offsets
   error_strings,		       %  error strings local to the function
};

private variable Error_Info_Map = Assoc_Type[Struct_Type];

private define lookup_error_info (func)
{
   if (assoc_key_exists (Error_Info_Map, func))
     return Error_Info_Map[func];
   return NULL;
}

private define error_handler ();
private define error_handler (func, code)
{
   variable info = lookup_error_info (func);
   if (info == NULL)
     throw ApplicationError, "No error info registered for $func (error=$code)"$;

   variable errstr;
   EXIT_BLOCK
     {
	return "$func:$errstr"$;
     }

   foreach (info.called)
     {
	variable called_error = ();
	if (code < 0)
	  {
	     if (called_error.offset > 0)
	       continue;
	     errstr = error_handler (called_error.func, -code);
	     return;
	  }

	if (code > called_error.offset)
	  {
	     errstr = error_handler (called_error.func, code - called_error.offset);
	     return;
	  }
     }

   try
     {
	errstr = info.error_strings[code-1];
     }
   catch IndexError: errstr = "Unknown error $code"$;
}

private define add_error (func, called_list, error_strings)
{
   variable info = @Error_Info_Type;
   variable called = {};
   _for (0, length(called_list)-1, 2)
     {
	variable i = ();
	list_append (called,
		     struct {func=called_list[i+1], offset=called_list[i]});
     }

   info.called = called;
   if (error_strings == NULL)
     error_strings = String_Type[0];
   info.error_strings = error_strings;

   Error_Info_Map[func] = info;
}

$2 = {20, "place", 10, "make_object",};
foreach $1 ({"app_star", "virtual_star", "astro_star"})
{
   add_error ($1, $2, NULL);
}

foreach $1 ({"topo_star", "local_star", "topo_planet", "local_planet"})
{
   add_error ($1, $2, ["Invalid value of 'where' in structure 'location'"]);
}

$2 = {10, "place"};
foreach $1 ({"app_planet", "virtual_planet", "astro_planet"})
{
   add_error ($1, $2, ["Invalid value of 'type' in structure 'ss_body'"]);
}


add_error("mean_star", {20, "app_star", 10, "vector2radec"},
	  ["Iterative process did not converge after 30 iterations",
	   "length of dummy star name out of bounds",
	   "length of dummy catalog name out of bounds",
	  ]);

add_error ("place",
	   {
	      90, "cio_basis", 80, "cio_location", 70, "grav_def", 50, "light_time",
		40, "geo_posvel", 10, "ephemeris",
	   },
	   ["invalid value of 'coord_sys'",
	    "invalid value of 'accuracy'",
	    "Earth is the observed object, and the observer is"
	    + " either at the geocenter or on the Earth's surface",
	   ]);

add_error ("equ2ecl", {}, ["invalid value of 'coord_sys'"]);
add_error ("equ2ecl_vec",  {}, ["invalid value of 'coord_sys'"]);
add_error ("ecl2equ_vec",  {}, ["invalid value of 'coord_sys'"]);

add_error ("gcrs2equ", {20, "cio_basis", 10, "cio_location", -1, "vector2radec"}, NULL);
add_error ("sidereal_time", {10, "cio_location",}, %cio_rai in docs!!!
	   ["invalid value of 'accuracy'", "invalid value of 'method'"]);
add_error ("ter2cel", {20, "cio_basis", 10, "cio_location",},
	   ["invalid value of 'accuracy'", "invalid value of 'method'"]);
add_error ("cel2ter", {20, "cio_basis", 10, "cio_location",},
	   ["invalid value of 'accuracy'", "invalid value of 'method'"]);

add_error ("cel_pole", {}, ["Invalid value of 'type'"]);
add_error ("geo_posvel", {}, ["invalid value of 'accuracy'"]);

add_error ("light_time", {10, "solarsystem"},
	   ["algorithm failed to converge after 10 iterations"]);
add_error ("grav_def", {30, "make_object", 0, "ephemeris"}, NULL);
add_error ("precession", {}, ["Precession not to or from J2000.0; 'jd_tdb1' or 'jd_tdb2' not 2451545.0"]);

add_error ("vector2radec", {},
	   ["All vector components are zero; 'ra' and 'dec' are indeterminate",
	    "Both pos[0] and pos[1] are zero, but pos[2] is nonzero; 'ra' is indeterminate",
	   ]);
add_error ("cio_ra", {20, "cio_basis", 10, "cio_location",},
	   ["invalid value of 'accuracy'"]);

add_error ("cio_location", {10, "cio_array"},
	   ["unable to allocate memory for the 'cio' array"]);

add_error ("cio_basis", {}, ["invalid value of input variable 'ref_sys'"]);

add_error ("cio_array", {},
	   ["error opening the 'cio_ra.bin' file.",
	    "'jd_tdb' not in the range of the CIO file.",
	    "n_pts' out of range",
	    "unable to allocate memory for the internal 't' array.",
	    "unable to allocate memory for the internal 'ra' array.",
	    "'jd_tdb' is too close to either end of the CIO file;"
	    + " unable to put 'n_pts' data points into the output structure.",
	   ]);
add_error ("ephemeris", {20, "readeph", 10, "solarsystem"},
	   ["Invalid value of 'origin'",
	    "Invalid value of 'type' in 'cel_obj'",
	    "Unable to allocate memory",
	   ]);
add_error ("transform_cat", {},
	   ["Invalid value of an input date for option 2 or 3 (see Note 1)",
	    "length of 'newcat_id' out of bounds.",
	   ]);
add_error ("make_object", {},
	   ["invalid value of 'type'",
	    "'number' out of range",
	    "Initialization of 'cel_obj' failed (object name)",
	    "Initialization of 'cel_obj' failed (catalog name).",
	    "'name' is out of string bounds.",
	   ]);

add_error ("readeph", {},
	   ["?1", "?2", "?3", "?4", "?5", "?6", "?7", "?8",
	    "This function is a stub for the real readeph and should not be called"];
	  );
add_error ("solarsystem", {},
	   ["illegal value for body",
	    "illegal value for origin",
	   ]);
$2 = "error reading from file header.";
add_error ("ephem_open", {},
	   [
	    "file does not exist/not found.",
	    $2, $2, $2, $2, $2, $2, $2, $2, $2,
	    "unable to set record length; ephemeris (DE number) not in look-up table.",
	   ]);

private define novas_error_callback (func, code)
{
   variable msg = error_handler (func, code);
   throw RunTimeError, msg;
}

novas_set_throw_callback (&novas_error_callback);

% The make_* functions are better implemented in slang.
define make_on_surface (lat, lon, height, temperature, pressure)
{
   return struct
     {
	latitude = lat,
	longitude = lon,
	height = height ,
	temperature = temperature,
	pressure = pressure,
     };
}

define make_in_space (sc_pos, sc_vel)
{
   return struct
     {
	sc_pos = sc_pos, sc_vel = sc_vel,
     };
}

define make_observer (where, obs_surface, obs_space)
{
   return struct
     {
	where = where,
	on_surf = obs_surface,
	near_earth = obs_space,
     };
}

define make_observer_on_surface (lat, lon, height, temperature, pressure)
{
   variable obs_surface = make_on_surface (lat, lon, height, temperature, pressure);
   variable obs_space = make_in_space ([0,0,0], [0,0,0]);
   return make_observer (1, obs_surface, obs_space);
}

define make_observer_in_space (pos, vel)
{
   variable obs_surface = make_on_surface (0, 0, 0, 0, 0);
   variable obs_space = make_in_space (pos, vel);
   return make_observer (2, obs_surface, obs_space);
}

define make_observer_at_geocenter ()
{
   variable obs_surface = make_on_surface (0, 0, 0, 0, 0);
   variable obs_space = make_in_space ([0,0,0],[0,0,0]);
   return make_observer (0, obs_surface, obs_space);
}

define make_cat_entry (starname, catalog, star_num, ra, dec,
		       pm_ra, pm_dec, parallax, rad_vel)
{
   return struct
     {
	starname = starname,
	catalog = catalog,
	starnumber = star_num,
	ra = ra,
	dec = dec,
	promora = pm_ra,
	promodec = pm_dec,
	parallax = parallax,
	radialvelocity = rad_vel,
     };
}

define make_object (type, number, name, star_data)
{
   return struct
     {
	type = type,
	number = number,
	name = name,
	star = star_data,
     };
}

