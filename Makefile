#####################################################################
#
#  Name:         Makefile
#  Created by:   Michael Mendenhall
#
#  Contents:     Makefile for RotationShield
#
#####################################################################
# export ROTSHIELD_VIS=1; make clean; make -j8; mv RotationShield RotationShieldVis; export ROTSHIELD_VIS=; make clean; make -j8

#############
# Some things you may need to change:
#############

# base path to usual unixy folders (include, lib, etc),
# probably just '/usr/'... unless you are using Fink on Mac OS X and have your stuff in '/sw/'
OS_DIR = /usr/

ifdef ROTSHIELD_MACOS
	OS_DIR = /sw/
	GL_PATH = /usr/X11
	GL_INCLUDES = -I/sw/include/ -I$(GL_PATH)/include/
	GL_LINKER_ARGS = $(GL_PATH)/lib/libGL.dylib /sw/lib/libglut.dylib
else
	GL_PATH = OS_DIR
	GL_INCLUDES = -I$(GL_PATH)/include/
	GL_LINKER_ARGS = -L$(GL_PATH)/lib/ -lGL -lglut
endif

# compiler command to use
CC = cc
CXX = g++

# include directories for GSL, FFTW, etc.
BASE_INCLUDE_DIRS  = -I$(OS_DIR)/include -IExtra_Headers
# lib dir flags for GSL, FFTW, etc.
BASE_LIB_DIRS  = -L$(OS_DIR)/lib

# optmization
GCC_OPTIMIZATION_LEVEL = 3

CPPFLAGS = -g -std=c++0x $(BUILDARCH) -O$(GCC_OPTIMIZATION_LEVEL) -pedantic -Wall -Wextra \
	-I. -IMathUtils -IGeometry -IFieldSource -ISolver -IBuilder -IStudies -IIO $(BASE_INCLUDE_DIRS)

#LDFLAGS += -lblas

ifdef ROTSHIELD_VIS
	CPPFLAGS += -DWITH_OPENGL $(GL_INCLUDES) 
	LDFLAGS += $(GL_LINKER_ARGS) -lpthread
endif

#############
# Everything below here "should" work without modification
#############

VPATH = ./:MathUtils/:Geometry/:FieldSource/:Solver/:Builder/:Studies/:IO/

# things to build
obj_IO = Visr.o strutils.o ControlMenu.o QFile.o SMExcept.o PathUtils.o VisSurface.o ProgressBar.o

obj_MathUtils = Geometry.o Integrator.o RefCounter.o linmin.o CMatrix.o InterpolationHelper.o BicubicGrid.o

obj_Geometry = Angles.o SurfaceGeometry.o SurfaceProfiles.o

obj_FieldSource = FieldSource.o MixedSource.o InfiniteLineSource.o LineSource.o InfinitePlaneSource.o \
	SurfaceSource.o SurfaceCurrentSource.o PlanarElement.o PlaneSource.o FieldEstimator2D.o

obj_Solver = ReactiveSet.o InterpolatingRS.o MagRS.o SurfaceCurrentRS.o GenericSolver.o SymmetricSolver.o

obj_Builder = CosThetaBuilder.o SurfacelCyl.o FieldAdaptiveSurface.o

obj_Studies = FieldAnalyzer.o tests.o Studies.o

objects = $(obj_IO) $(obj_MathUtils) $(obj_Geometry) $(obj_FieldSource) $(obj_Solver) $(obj_Builder) $(obj_Studies)
	
RotationShield : main.cpp $(objects)
	$(CXX) main.cpp $(objects) -o RotationShield $(CPPFLAGS) $(LDFLAGS) $(BASE_LIB_DIRS) -lgsl -lfftw3 -lgslcblas
	

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
	