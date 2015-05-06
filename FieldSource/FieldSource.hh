/* 
 * FieldSource.hh, part of the RotationShield program
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

/// \file FieldSource.hh \brief Base class for magnetic field sources

#ifndef FIELDSOURCE_HH
/// Make sure this header is only loaded once
#define FIELDSOURCE_HH

#include "Geometry.hh"
#include "SurfaceGeometry.hh"
#include "Typedefs.hh"
#include "RefCounter.hh"
#include "Integrator.hh"
#include "Matrix.hh"

/// Base (virtual) class for magnetic field sources due to current distributions
class FieldSource: public RefCounter, public Visualizable {
public:
    /// Constructor
    FieldSource(const std::string& nm = "FieldSource"): RefCounter(nm) {}
    /// Destructor
    virtual ~FieldSource() {}
    
    // Subclass me!
    //=====================================
    /// Magnetic field at a specified point
    virtual vec3 fieldAt(const vec3& v) const = 0;
    //=====================================
    
    /// Magnetic field at a point, with interaction matrix (optionally subclass for improved numerics)
    virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M2) const { return M2*fieldAt(v); }
    /// Magnetic field at a point, with interaction matrix (optionally subclass for improved numerics)
    virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M3) const { return M3*fieldAt(v); }
    /// Magnetic field averaged over specified line
    virtual vec3 fieldOverLine(Line l) const;
    /// Magnetic field averaged over specified plane
    virtual vec3 fieldOverPlane(Plane pl) const;
    /// Calculate effective net current of field source
    virtual vec3 net_current() const { return vec3(); }
    /// Print info to stdout
    virtual void display() const { printf("[FieldSource]\n"); }
    
    /// field averaged at or near surface
    virtual vec3 field_near(const SurfaceGeometry& S, double dh = 0, vec2 ll = vec2(0,0), vec2 ur = vec2(1,1)) const;
    /// RMS (component-wise) field strength averaged at or near surface
    virtual vec3 field_RMS_near(const SurfaceGeometry& S, double dh = 0, vec2 ll = vec2(0,0), vec2 ur = vec2(1,1)) const;
        
private:
    static Integrator lineIntegrator;   ///< Integrator for averaging over lines
    static Integrator planeIntegrator;  ///< Integrator for averaging over planes
};

#endif
