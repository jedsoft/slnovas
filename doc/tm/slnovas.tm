#% -*- mode: tm; mode: fold -*-

#%{{{Macros 

#i linuxdoc.tm
#d it#1 <it>$1</it>

#d slang \bf{S-lang}
#d exmp#1 \tt{$1}
#d var#1 \tt{$1}

#d ivar#1 \tt{$1}
#d ifun#1 \tt{$1}
#d cvar#1 \tt{$1}
#d cfun#1 \tt{$1}
#d svar#1 \tt{$1}
#d sfun#1 \tt{$1}
#d icon#1 \tt{$1}
#d dtype#1 \tt{$1}
#d exc#1 \tt{$1}

#d chapter#1 <chapt>$1<p>
#d preface <preface>
#d tag#1 <tag>$1</tag>

#d function#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d variable#1 \sect{<bf>$1</bf>\label{$1}}<descrip>
#d function_sect#1 \sect{$1}
#d begin_constant_sect#1 \sect{$1}<itemize>
#d constant#1 <item><tt>$1</tt>
#d end_constant_sect </itemize>

#d synopsis#1 <tag> Synopsis </tag> $1
#d keywords#1 <tag> Keywords </tag> $1
#d usage#1 <tag> Usage </tag> <tt>$1</tt>
#d description <tag> Description </tag>
#d qualifiers <tag> Qualifiers </tag>
#d example <tag> Example </tag>
#d notes <tag> Notes </tag>
#d seealso#1 <tag> See Also </tag> <tt>\linuxdoc_list_to_ref{$1}</tt>
#d done </descrip><p>
#d -1 <tt>-1</tt>
#d 0 <tt>0</tt>
#d 1 <tt>1</tt>
#d 2 <tt>2</tt>
#d 3 <tt>3</tt>
#d 4 <tt>4</tt>
#d 5 <tt>5</tt>
#d 6 <tt>6</tt>
#d 7 <tt>7</tt>
#d 8 <tt>8</tt>
#d 9 <tt>9</tt>
#d NULL <tt>NULL</tt>
#d documentstyle book

#%}}}

#d module#1 \tt{$1}
#d file#1 \tt{$1}
#d slang-documentation \
 \url{http://www.s-lang.org/doc/html/slang.html}{S-Lang documentation}

\linuxdoc
\begin{\documentstyle}

\title S-Lang NOVAS Module Reference
\author John E. Davis, \tt{jed@jedsoft.org}
\date \__today__

#i local.tm

\toc

\chapter{Introduction to the NOVAS Module}

Coordinate systems

GCRS:

The Geocentric Celestial Reference System (GCRS) is the system
appropriate for describing the rotation of the Earth, the orbits of
Earth satellites, and geodetic quantities such as instrument locations
and baselines. The directions of astronomical objects as seen from the
geocenter can also be expressed in the GCRS. The analysis of precise
observations inevitably involves quantities expressed in both systems
and the transformations between them.

The Barycentric Celestial Reference System (BCRS) is typically used for
the basic ephemerides of solar system objects and astrometric
reference data on galactic and extragalactic objects, i.e., the data
in astrometric star catalogs.

If the orientation of the BCRS axes in space is specified, the
orientation of the GCRS axes then follows from the relativistic
transformation between the two systems. The orientation of the BCRS is
given by what is called the International Celestial Reference System
(ICRS).

The ICRS is a triad of coordinate axes with their origin at the solar
system barycenter and with axis directions effectively defined by the
adopted coordinates of about 200 extragalactic radio sources observed
by VLBI (listed in Section H of The Astronomical Almanac). The
abbreviations BCRS and ICRS are often used interchangeably, because
the two concepts are so closely related: the ICRS is the orientation
of the BCRS while the BCRS is the metric for the ICRS.

The extragalactic radio sources that define the ICRS orientation are
assumed to have no observable intrinsic angular motions. Thus, the
ICRS is a ``space-fixed''system (more precisely, a
kinematically non-rotating system) and, as such, it has no associated
epoch--- its axes always point in the same directions with
respect to distant galaxies. However, the ICRS was set up to
approximate the conventional system defined by the Earth's
mean equator and equinox of epoch J2000.0; the alignment difference is
at the 0.02-arcsecond level, which is negligible for many
applications.

\chapter{NOVAS Module Function Reference}
#i slnovasfuns.tm

\end{\documentstyle}
