rotationshield
==============

magnetic field calculation inside rotationally symmetric magnetic shielding
by M. P. Mendenhall, 2007-2013

compile library dependencies: gsl fftw3 gslcblas

optional visualization dependencies: X11 libGL libglut
if you're feeling lucky, you can try compiling after setting
export ROTSHIELD_VIS=1

to compile:
make clean; make -j8

to run:
set your ROTSHIELD_OUT environment variable to indicate directory for output
run 'Rotationshield' executable produced by make for menu of options

example command lines:
bare coil:
./RotationShield coil geom 15 2.5 0.40 x meas svgrd 0 run foo x

shielded coil:
./RotationShield coil geom 15 2.5 0.40 x shield geom 2.9 0.47 grid 10 20 128 x meas svgrd 0 run foo x
