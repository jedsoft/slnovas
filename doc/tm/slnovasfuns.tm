#d ACCURACY \exmp{accuracy} specifies the accuracy of the calculation: 0 for full accuracy, 1 for reduced accuracy.
#d DELTAT  \exmp{delta_t} is the difference between in seconds between TT and UT1 at the specified Julian date

\function{app_star}
\synopsis{Compute the apparant place of a star}
\usage{(ra, dec) = app_star(jd_tt, cat_entry, accuracy)}
\description
 Computes the apparent place of a star at date \exmp{jd_tt}, given its
 catalog mean place, proper motion, parallax, and radial velocity.
 The second parameter is the catalog entry for the star, and \ACCURACY.

 The function returns the (ra, dec) of the object referred to the true
 equator and equinox of the specified date.  The right ascension is
 given in hours.
\seealso{virtual_star}
\done


\function{virtual_star}
\synopsis{Compute the virtual place of a star}
\usage{(ra, dec) = virtual_star(jd_tt, cat_entry, accuracy)}
\description
 Computes the virtual place of a star at date \exmp{jd_tt}, given its
 catalog mean place, proper motion, parallax, and radial velocity.
 The second parameter is the catalog entry for the star, and
 \exmp{accuracy} specifies the accuracy of the calculation: 0 for full
 accuracy, 1 for reduced accuracy.

 The function returns the (ra, dec) of the object referred to the GCRS.
 The right ascension is given in hour and the declination in degrees.
\seealso{app_star}
\done

\function{astro_star}
\synopsis{Compute the astrometric place of a star}
\usage{(ra, dec) = astro_star(jd_tt, cat_entry, accuracy)}
\description
 Computes the astrometric place of a star at date \exmp{jd_tt}, given its
 catalog mean place, proper motion, parallax, and radial velocity.
 The second parameter is the catalog entry for the star, and
 \exmp{accuracy} specifies the accuracy of the calculation: 0 for full
 accuracy, 1 for reduced accuracy.

 The function returns the (ra, dec) of the object referred to the
 ICRS, without light deflection or aberration. The right ascension is
 given in hours.
\seealso{app_star}
\done

\function{app_planet}
\synopsis{Compute the apparent place of a solar system body}
\usage{(ra, dec, dis) = app_planet (jd_tt, body, accuracy)}
\description
  This function computes the apparent place of the specified solar
  system body at the specified Julian date.  If \exmp{accuracy} is 0,
  the full accuracy calculation will be performed, and if 1, a reduced
  accuracy calculation will be done.

  The function returns the (ra, dec, dist) of the apparant position of the
  body with respect to the true equator and equinox of \exmp{jd_tt}.
  The right ascension is given in hours and the declination in degrees.
  The distance is in units of AU.
\seealso{app_star}
\done

\function{virtual_planet}
\synopsis{Compute the virtual place of a solar system body}
\usage{(ra, dec, dis) = virtual_planet (jd_tt, body, accuracy)}
\description
  This function computes the virtual place of the specified solar
  system body at the specified Julian date.  If \exmp{accuracy} is 0,
  the full accuracy calculation will be performed, and if 1, a reduced
  accuracy calculation will be done.

  The function returns the (ra, dec, dist) of the apparant position of
  the body with respect to the GCRS. The right ascension is given in
  hours and the declination in degrees.  The distance is in units of
  AU.
\seealso{app_star}
\done

\function{astro_planet}
\synopsis{Compute the astrometric place of a solar system body}
\usage{(ra, dec, dis) = astro_planet (jd_tt, body, accuracy)}
\description
  This function computes the astrometric place of the specified solar
  system body at the specified Julian date.  If \exmp{accuracy} is 0,
  the full accuracy calculation will be performed, and if 1, a reduced
  accuracy calculation will be done.

  The function returns the (ra, dec, dist) of the apparant position of the
  body with respect to the ICRS without light deflection or aberration
  The right ascension is given in hours and the declination in
  degrees.  The distance is in units of AU.
\notes
\seealso{}
\done


\function{topo_star}
\synopsis{Compute the topocentric place of a star}
\usage{(ra, dec) = topo_star(jd_tt, delta_t, star, position, accuracy)}
\description
 Computes the topocentric place of a star at Julian date \exmp{jd_tt},
 given its catalog mean place, proper motion, parallax, and radial
 velocity. \DELTAT, \exmp{star} is the
 catalog entry for the object in the ICRS, and \exmp{position} is the
 surface object for the position of the observer.  If \exmp{accuracy}
 is 0, a full accuracy calculation will be performed, otherwise if 1,
 a reduced accuracy calculation will be done.

 The function returns the (ra in hours, dec in degrees) with respect
 to the true equator and equinox of date \exmp{jd_tt}.
\seealso{local_star}
\done


\function{local_star}
\synopsis{Compute the local place of a star}
\usage{(ra, dec) = local_star(jd_tt, delta_t, star, position, accuracy)}
\description
 Computes the local place of a star at Julian date \exmp{jd_tt},
 given its catalog mean place, proper motion, parallax, and radial
 velocity.  \exmp{delta_t} is the difference between in seconds
 between TT and UT1 at the specified Julian date, \exmp{star} is the
 catalog entry for the object in the ICRS, and \exmp{position} is the
 surface object for the position of the observer.  If \exmp{accuracy}
 is 0, a full accuracy calculation will be performed, otherwise if 1,
 a reduced accuracy calculation will be done.

 The function returns the local (ra in hours, dec in degrees) with respect
 to the local GCRS.
\seealso{topo_star}
\done


\function{topo_planet}
\synopsis{Compute the local place of a solar system body}
\usage{(ra, dec, dis) = topo_planet(jd_tt, ss_body, delta_t, position, accuracy)}
\description
 Computes the topocentric place of a solar system body at the Julian
 date \exmp{jd_tt}.  \exmp{delta_t} is the difference
 between in seconds between TT and UT1 at the specified Julian date,
 \exmp{ss_body} specifies the solar system object, and
 \exmp{position} is the surface object for the position of the
 observer.  If \exmp{accuracy} is 0, a full accuracy calculation will
 be performed, otherwise if 1, a reduced accuracy calculation will be
 done.

 The function returns the local (ra in hours, dec in degrees, dist in AU) with respect
 to the true equator and equinox of the date.
\seealso{local_planet, astro_planet, topo_star}
\done

\function{local_planet}
\synopsis{Compute the local place of a solar system body}
\usage{(ra, dec, dis) = local_planet(jd_tt, ss_body, delta_t, position, accuracy)}
\description
 Computes the local place of a solar system body at the Julian
 date \exmp{jd_tt}.  \exmp{delta_t} is the difference
 between in seconds between TT and UT1 at the specified Julian date,
 \exmp{ss_body} specifies the solar system object, and
 \exmp{position} is the surface object for the position of the
 observer.  If \exmp{accuracy} is 0, a full accuracy calculation will
 be performed, otherwise if 1, a reduced accuracy calculation will be
 done.

 The function returns the local (ra in hours, dec in degrees, dist in
 AU) with respect to the local GCRS.
\seealso{topo_planet, virtual_planet, local_star}
\done

\function{mean_star}
\synopsis{}
\usage{(ira, idec) = mean_star(jd_tt, ra, dec, accuracy)}
\description
  Computes the ICRS position of a star, given its apparent place at
  date \exmp{jd_tt}.  Proper motion, parallax and radial velocity are
  assumed to be zero.  The input arguments \exmp{(ra, dec)} are the
  apparant right ascension in hours and declination in degrees with
  respect to the true equator and equinox of the date. If
  \exmp{accuracy} is 0, a full accuracy calculation will be performed,
  otherwise if 1, a reduced accuracy calculation will be done.

  The function returns the ICRS (ra in hours, dec in degrees).
\example
\notes
\seealso{}
\done

\function{place}
\synopsis{Compute the apparent direction of a star or solar system body}
\usage{struct sky_pos = place(jd_tt, cel_object, location, delta_t, coord_sys, accuracy)}
\description
  This function computes the apparent direction of a star or solar
  system body at the Julian date \exmp{jd_tt} and in specified coordinate
  system.  The \exmp{cel_object} parameter specifies the celestial
  object, the \exmp{location} parameter specifies the observers
  location, \exmp{delta_t} is the difference between TT and UT1 and
  \exmp{jd_tt}, and \exmp{coord_sys} is one of the following:
#v+
    0 ... GCRS or "local GCRS"
    1 ... true equator and equinox of date
    2 ... true equator and CIO of date
    3 ... astrometric coordinates, i.e., without light deflection or aberration.
#v-
  If \exmp{accuracy} is 0, a full accuracy calculation will be
  performed, otherwise a reduced accuracy one will be done.

  The function returns a structure giving the position of the object
  in the sky with respect to the specified coordinate system.
\notes
   Values of \exmp{location.where} and \exmp{coord_sys} dictate the various
   standard kinds of place:
#v+
   location.where = 0 and coord_sys = 1: apparent place
   location.where = 1 and coord_sys = 1: topocentric place
   location.where = 0 and coord_sys = 0: virtual place
   location.where = 1 and coord_sys = 0: local place
   location.where = 0 and coord_sys = 3: astrometric place
   location.where = 1 and coord_sys = 3: topocentric astrometric place
#v-

 The value of \exmp{delta_t} is used only when \exmp{location.where}
 equals 1 or 2 (observer is on surface of Earth or in a
 near-Earth satellite).
\seealso{topo_planet}
\done

\function{julian_date}
\synopsis{}
\usage{jd = julian_date(year, month, day, hour)}
\description
\example
\notes
\seealso{}
\done


\function{equ2gal}
\synopsis{Convert ICRS ra/dec to galactic lon/lat}
\usage{(ra, dec) = equ2gal(rai, deci)}
\description
  This function converts the specified ICRS right ascension and
  declination to galactic longitude and latitude and returns the result.
\seealso{equ2ecl, equ2hor}
\done

\function{equ2ecl}
\synopsis{Convert ra/dec to ecliptic lon/lat}
\usage{(elon, elat) = equ2ecl(jd_tt, coord_sys, accuracy, ra, dec)}
\description
 This function converts right ascension and declination to ecliptic longitude
 and latitude.  The \exmp{jd_tt} argument gives the Julian date of the
 equator, equinox, and ecliptic used.  The \exmp{coord_sys} parameter
 specifies the coordinate system:
#v+
  0 ... mean equator and equinox of date 'jd_tt'
  1 ... true equator and equinox of date 'jd_tt'
  2 ... ICRS
  (ecliptic is always the mean plane)
#v-
  If \exmp{accuracy} is 0, a full accuracy calculation will be
  performed, otherwise a reduced accuracy one will be done.
\notes
\seealso{equ2gal, equ2hor, equ2ecl_vec}
\done

\function{equ2hor}
\synopsis{transforms topocentric ra/dec to zenith distance and azimuth}
\usage{(zd,az,rar,decr)=equ2hor(jd_ut1, delta_t, accuracy, xp, yp, location, ra, dec, ref_option)}
\description
  This function transforms topocentric right ascension and
  declination to zenith distance and azimuth.  It uses a method
  that properly accounts for polar motion, which is significant at
  the sub-arcsecond level.  This function can also adjust
  coordinates for atmospheric refraction.

  \exmp{jd_ut1} specifies the UT1 Julian date in seconds, and \DELTAT.
  The parameters \exmp{(xp,yp)} give the conventionally defined x and
  y coordinates of the intermediate pole with respect to the ITRS
  reference pole in arc-seconds.  The observer's surface location is
  given by the \exmp{location} argument, and \exmp{(ra,dec)} specifies
  the topocentric right ascension and declination of the object of
  interest with respect to the true equator and equinox of the date.
  The \exmp{ref_option} parameter specifies how refraction should be
  handled:
#v+
    0 ... no refraction
    1 ... include refraction, using 'standard' atmospheric conditions.
    2 ... include refraction, using atmospheric parameters input in
      the 'location' structure.
#v-

  The function returns the following 4 values:
#v+
  zd: Topocentric zenith distance in degrees, affected by
       refraction if 'ref_option' is non-zero.
  az: Topocentric azimuth (measured east from north) in degrees.
 rar: Topocentric right ascension of object of interest, in hours,
       referred to true equator and equinox of date, affected by
       refraction if 'ref_option' is non-zero.
decr: Topocentric declination of object of interest, in degrees,
       referred to true equator and equinox of date, affected by
       refraction if 'ref_option' is non-zero.
#v-
\seealso{equ2gal}
\done

\function{sidereal_time}
\synopsis{Compute the Greenwich mean or apparent sidereal time}
\usage{t = sidereal_time(jd_high, jd_low, delta_t, gst_type, method, accuracy)}
\description
  This function computes the Greenwich sidereal time, either mean or
  apparent, at Julian date date \exmp{jd_high+jd_low}.
  The parameter \DELTAT.
  If \exmp{gst_type} is 0, the Greenwich mean sidereal time will be
  computed, and if \exmp{gst_type} is 1, the apparent time will be
  computed.  If \exmp{method} is 0, a CIO-based method will be used,
  otherwise a equinox-based method will be used.  \ACCURACY.

  The function returns the Greenwich (mean or apparent) sidereal time
  in hours.
\seealso{era, julian_date}
\done


\function{era}
\synopsis{Compute the value of the Earth Rotation Angle}
\usage{theta = era(jd_high, jd_low)}
\description
  This function returns the value of the Earth Rotation Angle (theta)
  for a given UT1 Julian date.  The expression used is taken from the
  note to IAU Resolution B1.8 of 2000.  The parameters
  \exmp{jd_high/jd_low} specify the high/low order parts of the UT1
  Julian date.

  The function returns the Earth Rotation Angle in degrees.
\seealso{julian_date}
\done


\function{ephemeris}
\synopsis{}
\usage{ephemeris()}
\description
\example
\notes
\seealso{}
\done


\function{equ2ecl_vec}
\synopsis{Convert an equatorial position vector to an ecliptic position vector}
\usage{v2 = equ2ecl_vec(jd_tt, coord_sys, accuracy, v1)}
\description
 This function Converts an equatorial position vector \exmp{v1} to an
 ecliptic position vector.  The \exmp{jd_tt} parameter gives the
 Julian date of the equator, equinox, and ecliptic used for the
 coordinates, and coord_sys specifies the coordinate system to use:
#v+
  0 ... mean equator and equinox of date 'jd_tt'
  1 ... true equator and equinox of date 'jd_tt'
  2 ... ICRS
  (ecliptic is always the mean plane)
#v-
 \ACCURACY.

 The function returns the transformed vector referred to specified
 ecliptic and equinox of date.
\notes
  To convert an ICRS vector to an ecliptic vector (mean ecliptic
  and equinox of J2000.0 only), set \exmp{coord_sys = 2}; the value
  of \exmp{jd_tt} can be set to anything, since J2000.0 is assumed.

  Except for the input to this case, all vectors are assumed
  to be with respect to a dynamical system.
\seealso{ecl2equ_vec}
\done


\function{ecl2equ_vec}
\synopsis{}
\usage{ecl2equ_vec()}
\description
\example
\notes
\seealso{}
\done


\function{gcrs2equ}
\synopsis{}
\usage{gcrs2equ()}
\description
\example
\notes
\seealso{}
\done


\function{vector2radec}
\synopsis{}
\usage{vector2radec()}
\description
\example
\notes
\seealso{}
\done


\function{ter2cel}
\synopsis{Rotate a vector from the terrestrial to the celestial system}
\usage{v_cel = ter2cel(jd_ut_high, jd_ut_low, delta_t, method, accuracy, option, xp, yp, v_ter)}
\description
  This function rotates a vector from the terrestrial to the
  celestial system.  Specifically, it transforms a vector in the
  ITRS (rotating earth-fixed system) to the GCRS (a local space-fixed
  system) by applying rotations for polar motion, Earth
  rotation, nutation, precession, and the dynamical-to-GCRS
  frame tie.

  The parameters \exmp{(jd_ut_high, jd_ut_low)} specify the high and
  low order parts of the UT1 Julian date.  The parameter \DELTAT.  If
  \exmp{method} is 0, a CIO-based method will be used, otherwise an
  equinox-based on will be used.  \ACCURACY.  If \exmp{option} is 0,
  the output will be referred to the GCRS axes, and if 1, it will be
  with respect to the equator and equinox of date.  The \exmp{(xp,yp)}
  parameters give the conventionally-defined (x,y) coordinates of the
  ITRS pole in arc-seconds.  The vector \exmp{v_ter} is the vector to be
  rotated, in geocentric equatorial rectangular coordinates,
  referred to ITRS axes (terrestrial system) in the normal case
  where 'option' = 0.

  The function returns a vector in geocentric equatorial rectangular
  coordinates, referred to GCRS axes (celestial system) or with
  respect to the equator and equinox of date, depending on
  \exmp{option}.
\notes
  \exmp{xp = yp = 0} means no polar motion transformation.

  The \exmp{option} flag only works for the equinox-based method.
\seealso{cel2ter}
\done

\function{cel2ter}
\synopsis{Rotate a vector from the celestial to the terrestrial system}
\usage{v_ter = cel2ter(jd_ut_high, jd_ut_low, delta_t, method, accuracy, option, xp, yp, v_cel)}
\description
  This function is the reverse of \ifun{ter2cel}.  See its
  documentation for a description of the input arguments.
\seealso{ter2cel}
\done

\function{ephem_open}
\synopsis{}
\usage{ephem_open()}
\description
\example
\notes
\seealso{}
\done


\function{ephem_close}
\synopsis{}
\usage{ephem_close()}
\description
\example
\notes
\seealso{}
\done


\function{aberration}
\synopsis{}
\usage{aberration()}
\description
\example
\notes
\seealso{}
\done


\function{bary2obs}
\synopsis{}
\usage{bary2obs()}
\description
\example
\notes
\seealso{}
\done


\function{cal_date}
\synopsis{}
\usage{cal_date()}
\description
\example
\notes
\seealso{}
\done


\function{cel_pole}
\synopsis{}
\usage{cel_pole()}
\description
\example
\notes
\seealso{}
\done


\function{cio_array}
\synopsis{}
\usage{cio_array()}
\description
\example
\notes
\seealso{}
\done


\function{cio_basis}
\synopsis{}
\usage{cio_basis()}
\description
\example
\notes
\seealso{}
\done


\function{cio_location}
\synopsis{}
\usage{cio_location()}
\description
\example
\notes
\seealso{}
\done


\function{cio_ra}
\synopsis{}
\usage{cio_ra()}
\description
\example
\notes
\seealso{}
\done


\function{d_light}
\synopsis{}
\usage{d_light()}
\description
\example
\notes
\seealso{}
\done


\function{ee_ct}
\synopsis{}
\usage{ee_ct()}
\description
\example
\notes
\seealso{}
\done


\function{e_tilt}
\synopsis{}
\usage{e_tilt()}
\description
\example
\notes
\seealso{}
\done


\function{frame_tie}
\synopsis{}
\usage{frame_tie()}
\description
\example
\notes
\seealso{}
\done


\function{fund_args}
\synopsis{}
\usage{fund_args()}
\description
\example
\notes
\seealso{}
\done


\function{geo_posvel}
\synopsis{}
\usage{geo_posvel()}
\description
\example
\notes
\seealso{}
\done


\function{grav_def}
\synopsis{}
\usage{grav_def()}
\description
\example
\notes
\seealso{}
\done


\function{grav_vec}
\synopsis{}
\usage{grav_vec()}
\description
\example
\notes
\seealso{}
\done


\function{ira_equinox}
\synopsis{}
\usage{ira_equinox()}
\description
\example
\notes
\seealso{}
\done


\function{light_time}
\synopsis{}
\usage{light_time()}
\description
\example
\notes
\seealso{}
\done


\function{limb_angle}
\synopsis{}
\usage{limb_angle()}
\description
\example
\notes
\seealso{}
\done


\function{mean_obliq}
\synopsis{}
\usage{mean_obliq()}
\description
\example
\notes
\seealso{}
\done


\function{norm_ang}
\synopsis{}
\usage{norm_ang()}
\description
\example
\notes
\seealso{}
\done


\function{nutation}
\synopsis{}
\usage{nutation()}
\description
\example
\notes
\seealso{}
\done


\function{nutation_angles}
\synopsis{}
\usage{nutation_angles()}
\description
\example
\notes
\seealso{}
\done


\function{planet_ephemeris}
\synopsis{}
\usage{planet_ephemeris()}
\description
\example
\notes
\seealso{}
\done


\function{precession}
\synopsis{}
\usage{precession()}
\description
\example
\notes
\seealso{}
\done


\function{proper_motion}
\synopsis{}
\usage{proper_motion()}
\description
\example
\notes
\seealso{}
\done


\function{radec2vector}
\synopsis{}
\usage{radec2vector()}
\description
\example
\notes
\seealso{}
\done


\function{rad_vel}
\synopsis{}
\usage{rad_vel()}
\description
\example
\notes
\seealso{}
\done


\function{refract}
\synopsis{}
\usage{refract()}
\description
\example
\notes
\seealso{}
\done


\function{spin}
\synopsis{Rotate a vector about the z-axis}
\usage{v2 = spin(angle, v1)}
\description
  This function rotates the vector \exmp{v1} about the z-axis by an
  \exmp{angle} degrees.  The rotated angle is returned.
\notes
  This corresponds to a right-handed rotation about the +z-axis.
\seealso{wobble}
\done


\function{starvectors}
\synopsis{}
\usage{starvectors()}
\description
\example
\notes
\seealso{}
\done


\function{tdb2tt}
\synopsis{}
\usage{tdb2tt()}
\description
\example
\notes
\seealso{}
\done


\function{terra}
\synopsis{}
\usage{terra()}
\description
\example
\notes
\seealso{}
\done


\function{transform_cat}
\synopsis{}
\usage{transform_cat()}
\description
\example
\notes
\seealso{}
\done


\function{transform_hip}
\synopsis{}
\usage{transform_hip()}
\description
\example
\notes
\seealso{}
\done


\function{wobble}
\synopsis{}
\usage{wobble()}
\description
\example
\notes
\seealso{}
\done


\function{novas_set_throw_callback}
\synopsis{}
\usage{novas_set_throw_callback()}
\description
\example
\notes
\seealso{}
\done
