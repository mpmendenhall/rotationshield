#####################################################################
#
#  Name:         Makefile
#  Created by:   Michael Mendenhall
#
#  Contents:     Makefile for RotationShield
#
#####################################################################

#############
# Some things you may need to change:
#############

# compiler command to use
CC = cc
CXX = g++

# base path to usual unixy folders (include, lib, etc),
# probably just '/'... unless you are using Fink on Mac OS X and have your stuff in '/sw/'
OS_DIR = /sw
# include directories for GSL, CLHEP, FFTW, etc.
BASE_INCLUDE_DIRS  = -I$(OS_DIR)/include
# lib dir flags for GSL, FFTW, etc.
BASE_LIB_DIRS  = -L$(OS_DIR)/lib
# architecture flag... need i386 for old 32-bit Fink depends
BUILDARCH = -arch i386
# locations of X11 and OpenGL-related stuff
X11_PATH = /Developer/SDKs/MacOSX10.6.sdk/usr/X11
GL_INCLUDES = -I$(X11_PATH)/include/
GLUT_LIB = $(X11_PATH)/lib/
# linker arguments needed for OpenGL/GLUT visualization window
GL_LINKER_ARGS = -L$(GLUT_LIB) -lGL -lglut

CPPFLAGS = -g $(BUILDARCH) -O$(GCC_OPTIMIZATION_LEVEL) -Wall -Wuninitialized $(BASE_INCLUDE_DIRS) $(BASE_LIB_DIRS) $(USERFLAGS) $(GL_INCLUDES) \
	-I. -IMathUtils -IFieldSource -ISolver -IBuilder -IStudies -IIO

#############
# Everything below here "should" work without modification
#############

VPATH = ./:MathUtils/:FieldSource/:Solver/:Builder/:Studies/:IO/

# things to build
obj_IO = Visr.o

obj_MathUtils = Geometry.o Integrator.o MiscUtils.o RefCounter.o

obj_FieldSource = FieldSource.o MixedSource.o InfiniteLineSource.o LineSource.o InfinitePlaneSource.o PlaneSource.o Boxel.o

obj_Solver = GenericSolver.o InteractionSolver.o SymmetricSolver.o

obj_Builder = ShieldBuilder.o

objects = $(obj_IO) $(obj_MathUtils) $(obj_FieldSource) $(obj_Solver) $(obj_Builder)
	
main :
	make RotationShield

RotationShield : main.cpp $(objects)
	$(CXX) main.cpp $(objects) -o RotationShield $(CPPFLAGS) $(GL_LINKER_ARGS) -lCLHEP -lgsl -lfftw3
	

#
# documentation via Doxygen
#
doc : latex/refman.pdf

latex/refman.pdf: latex/ 
	cd latex; make
latex/ : Doxyfile $(objects)
	doxygen

#
# cleanup
#
clean:
	-rm -rf *.o
	-rm -rf latex/
	-rm -rf html/
	-rm RotationShield
	