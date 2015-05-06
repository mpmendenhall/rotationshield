#####################################################################
#
# Makefile, part of the RotationShield program
# Copyright (c) 2007-2014 Michael P. Mendenhall
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
######################################################################

# export ROTSHIELD_VIS=1; make clean; make -j8; mv RotationShield RotationShield_Vis; export ROTSHIELD_VIS=; make clean; make -j8



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
	GL_PATH = $(OS_DIR)
	GL_INCLUDES = -I$(GL_PATH)/include/
	GL_LINKER_ARGS = -L$(GL_PATH)/lib/ -lGL -lglut
endif

# compiler command to use
CC = gcc
CXX = g++

# include directories for GSL, FFTW, etc.
BASE_INCLUDE_DIRS  = -I$(OS_DIR)/include -IExtra_Headers
# lib dir flags for GSL, FFTW, etc.
BASE_LIB_DIRS  = -L$(OS_DIR)/lib

# optmization
GCC_OPTIMIZATION_LEVEL = 3

CFLAGS = -g $(BUILDARCH) -O$(GCC_OPTIMIZATION_LEVEL) -pedantic -Wall -Wextra -Icubature-1.0 $(BASE_INCLUDE_DIRS)

CXXFLAGS = -g -std=c++0x $(BUILDARCH) -O$(GCC_OPTIMIZATION_LEVEL) -pedantic -Wall -Wextra \
	-I. -IMathUtils -IGeometry -IFieldSource -ISolver -IBuilder -IStudies -IIO -Icubature-1.0 $(BASE_INCLUDE_DIRS)

ifdef ROTSHIELD_LAPACKE
	# BLAS and LAPACK(E) libraries for matrix manipulation; also, need gfortran for LAPACK linking
	LDFLAGS += -llapacke -llapack -lblas -lgfortran
	BASE_LIB_DIRS += -L$(OS_DIR)/lib/lapack/ -L$(OS_DIR)/lib/gcc4.9/lib/
	CXXFLAGS += -DWITH_LAPACKE
endif

ifdef ROTSHIELD_VIS
	CXXFLAGS += -DWITH_OPENGL $(GL_INCLUDES) 
	LDFLAGS += $(GL_LINKER_ARGS) -lpthread
endif

#############
# Everything below here "should" work without modification
#############

VPATH = ./:MathUtils/:Geometry/:FieldSource/:Solver/:Builder/:Studies/:IO/:cubature-1.0

LDFLAGS += $(BASE_LIB_DIRS) -lgsl -lfftw3 -lgslcblas

# things to build
obj_IO = Visr.o strutils.o ControlMenu.o QFile.o SMExcept.o PathUtils.o VisSurface.o ProgressBar.o

obj_MathUtils = Geometry.o Integrator.o CMatrix.o LAPACKE_Matrix.o BlockCMat.o RefCounter.o \
	linmin.o InterpolationHelper.o BicubicGrid.o hcubature.o pcubature.o

obj_Geometry = Angles.o SurfaceGeometry.o SurfaceProfiles.o

obj_FieldSource = FieldSource.o MixedSource.o InfiniteLineSource.o LineSource.o InfinitePlaneSource.o \
	SurfaceSource.o SurfaceCurrentSource.o PlanarElement.o PlaneSource.o SymmetrizedSource.o DipoleSource.o FieldEstimator2D.o

obj_Solver = ReactiveSet.o InterpolatingRS.o MagRS.o SurfaceCurrentRS.o HoleDipolePerturbation.o \
	InteractionSolver.o GenericSolver.o SymmetricSolver.o MultiQuilibrator.o

obj_Builder = CosThetaBuilder.o SurfacelCyl.o FieldAdaptiveSurface.o

obj_Studies = FieldAnalyzer.o tests.o Studies.o

objects = $(obj_IO) $(obj_MathUtils) $(obj_Geometry) $(obj_FieldSource) $(obj_Solver) $(obj_Builder) $(obj_Studies)

all: RotationShield

libRotationShield.a: $(objects)
	ar rs libRotationShield.a $(objects)

RotationShield: main.cpp libRotationShield.a
	$(CXX) $(CXXFLAGS) main.cpp -L. -lRotationShield $(LDFLAGS) -o RotationShield


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
	-rm -f RotationShield
	-rm -f libRotationShield.a
	
