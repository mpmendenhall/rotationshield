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

shielded with superconducting endcaps:
./RotationShield coil geom 15 2.5 0.40 x shield geom 2.9 0.47 grid 10 20 128 ecap 6 0.25 0 both x meas svgrd 0 run foo x

half-scale model with single sided superconducting endcap:
./RotationShield coil geom 15 2.146 0.324 x   shield geom 2.159 0.362 grid 10 20 128    ecap 10 0.102 0 neg  mcap 0.051 0.051 neg x    cell range -.18 -.18 -0.4 .18 .18 0.4 grid 7 7 13 x    meas svgrd 1 run halfscale_sc_end x

same, without endcap:
./RotationShield coil geom 15 2.146 0.324 x   shield geom 2.159 0.362 grid 10 20 128 x   cell range -.18 -.18 -0.4 .18 .18 0.4 grid 7 7 13 x    meas svgrd 1 run halfscale_open_end x
