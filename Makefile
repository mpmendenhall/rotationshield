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

# base path to usual unixy folders (include, lib, etc),
# probably just '/usr/'... unless you are using Fink on Mac OS X and have your stuff in '/sw/'
OS_DIR = /usr/

# compiler command to use
CC = cc
CXX = g++

# include directories for GSL, CLHEP, FFTW, etc.
BASE_INCLUDE_DIRS  = -I$(OS_DIR)/include
# lib dir flags for GSL, FFTW, etc.
BASE_LIB_DIRS  = -L$(OS_DIR)/lib
# locations of X11 and OpenGL-related stuff
X11_PATH = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/X11
GL_INCLUDES = -I$(X11_PATH)/include/
GLUT_LIB = $(X11_PATH)/lib/
# linker arguments needed for OpenGL/GLUT visualization window
GL_LINKER_ARGS = -L$(GLUT_LIB) -lGL -lglut
# optmization
GCC_OPTIMIZATION_LEVEL = 3

CPPFLAGS = -g $(BUILDARCH) -O$(GCC_OPTIMIZATION_LEVEL)\
		-Wall -Wuninitialized -I. -IMathUtils -IFieldSource -ISolver -IBuilder -IStudies -IIO\
		$(BASE_INCLUDE_DIRS) $(USERFLAGS)


#############
# Everything below here "should" work without modification
#############

VPATH = ./:MathUtils/:FieldSource/:Solver/:Builder/:Studies/:IO/

# things to build
obj_IO = Visr.o

obj_MathUtils = Geometry.o Integrator.o MiscUtils.o RefCounter.o analysis.o linmin.o

obj_FieldSource = FieldSource.o MixedSource.o InfiniteLineSource.o LineSource.o InfinitePlaneSource.o PlaneSource.o Boxel.o

obj_Solver = GenericSolver.o InteractionSolver.o SymmetricSolver.o

obj_Builder = ShieldBuilder.o CosThetaBuilder.o

obj_Studies = tests.o

objects = $(obj_IO) $(obj_MathUtils) $(obj_FieldSource) $(obj_Solver) $(obj_Builder) $(obj_Studies)
	
RotationShield : main.cpp $(objects)
	$(CXX) main.cpp $(objects) -o RotationShield $(CPPFLAGS) $(BASE_LIB_DIRS) -lCLHEP -lgsl -lfftw3 -lgslcblas
	

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
	