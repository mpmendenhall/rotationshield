==============
RotationShield
==============

Linear matrix solver for systems with cyclic symmetry,
with application to magnetic fields in rotationally symmetric geometries.

=========================================================================

Copyright (c) 2007-2014 Michael P. Mendenhall

and including a re-distribution of "cubature," (under GPL v2 or later)
Copyright (c) 2005-2013 Steven G. Johnson; http://ab-initio.mit.edu/cubature/

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=========================================================================


compile library dependencies: gsl fftw3 gslcblas

optional visualization dependencies: X11 libGL libglut
if you're feeling lucky, you can try compiling after setting
export ROTSHIELD_VIS=1

optional SVD decomposition for symmetric solver; needs LAPACK, LAPACKE, BLAS, gfortran libraries available
export ROTSHIELD_LAPACKE=1
(and edit makefile as needed to find correct libraries)

to compile:
make clean; make -j8

to run:
set your ROTSHIELD_OUT environment variable to indicate directory for output

run 'Rotationshield' executable produced by make for menu of options.
	Note, you may abbreviate menu commands to the shortest unambiguous initial substring.

'Rotationshield' may also be scripted from the command line,
	provided with the same argument sequence used for menu control.
