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
run 'Rotationshield' executable produced by make for menu of options

----------------------
example command lines:
----------------------

# bare coil:
./RotationShield coil geom 15 2.5 0.40 x meas svgrd 0 run foo x

# shielded coil:
./RotationShield coil geom 15 2.5 0.40 x  shield geom 2.9 0.47 128  add 10000 7 15 nrd prd x  meas svgrd 0 run foo x

# shielded with superconducting endcaps:
./RotationShield coil geom 15 2.5 0.40 x  shield geom 2.9 0.47 128  add 10000 7 15 nrd prd  add 0 12 6 nax nrd  add 0 12 6 prd pax x  meas svgrd 0 run foo x

# half-scale model with single sided annular superconducting endcap:
./RotationShield coil geom 15 2.146 0.324 x  shield geom 2.159 0.362 128  add 10000 10 20 nrd prd  add 0 12 6 nax nrd  mod 0 0.102 0 0 x   cell range -.18 -.18 -0.4 .18 .18 0.4 grid 7 7 13 x    meas svgrd 1 run halfscale_sc_end x
