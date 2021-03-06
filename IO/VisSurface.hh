/* 
 * VisSurface.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/// \file VisSurface.hh \brief Visualize parameters mapped onto a surface
#ifndef VISSURFACE_HH
/// Make sure this header is only loaded once
#define VISSURFACE_HH

#include "Visr.hh"
#include "SurfaceGeometry.hh"
#include "BicubicGrid.hh"
#include "Color.hh"
#include "Typedefs.hh"

/// general purpose visualization surface
class VisSurface: public Visualizable {
public:
    /// constructor
    VisSurface(SurfaceGeometry* SG): mySurface(SG), vis_nx(100), vis_ny(100) { }
    /// destructor
    virtual ~VisSurface() { }
    
    /// visualization routine
    virtual void _visualize() const;
    
    SurfaceGeometry* mySurface;    ///< surface over which current is distributed
    unsigned int vis_nx;        ///< fineness of x grid
    unsigned int vis_ny;        ///< fineness oy y grid

    double (*f_height)(vec2, const void*);    ///< color generator
    void* hparams;
    vec4 (*f_color)(vec2, const void*);        ///< surface normal offset generator
    void* cparams;
};

/// surface for sampling+interpolating, to provide visualization data
class SamplerSurface {
public:
    /// constructor
    SamplerSurface() {}
    /// destructor
    virtual ~SamplerSurface() { clear_data(); }
    
    /// fill grid data from vector function
    void sample_func(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* fparams);
    /// evaluate sampled data
    mvec eval(vec2) const;
    /// evaluate one axis
    double eval(unsigned int n, vec2 l) const { assert(n<G.size()); return (*G[n])(l[0],l[1]); }
    
    /// delete previous grids
    void clear_data();
    /// set up data grids
    void make_grids(unsigned int nx, unsigned int ny, unsigned int ng);
    
    vector<BicubicGrid*> G;
    unsigned int NX;
    unsigned int NY;
    
};

/// utility function for height from one variable
double h_fromSampler1(vec2 l, const void* smplr);
/// utility function for height from vector magnitude
double h_fromSamplerMag(vec2 l, const void* smplr);
/// utility function for color from 2 variable direction
vec4 c_fromSampler2(vec2 l, const void* smplr);

/// visualize data passed to 2D integrator
void vis_integrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* params = NULL);

#endif
