rotationshield
==============

magnetic field calculation inside rotationally symmetric magnetic shielding
by M. P. Mendenhall, 2007-2014

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

----------------------
example command lines: TODO
----------------------
