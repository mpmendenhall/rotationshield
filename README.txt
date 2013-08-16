rotationshield
==============

magnetic field calculation inside rotationally symmetric magnetic shielding
by M. P. Mendenhall, 2007-2013

compile library dependencies: CLHEP gsl fftw3 gslcblas

optional visualization dependencies: X11 libGL libglut
if you're feeling lucky, you can try compiling after setting
export ROTSHIELD_VIS="1"

to compile:
make clean; make -j8

to run:
set your ROTSHIELD_OUT to indicate directory for output
run 'Rotationshield' executable produced by make for menu of options
