/* 
 * UniformField.hh, part of the RotationShield program
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

#ifndef UNIFORMFIELD_HH
/// Make sure this header is only loaded once
#define UNIFORMFIELD_HH

#include "FieldSource.hh"
#include <stdio.h>

/// A uniform magnetic field
class UniformField: public FieldSource {
public:
    /// Constructor
    UniformField(const vec3& f = vec3()): FieldSource("UniformField"), B(f) { mySymmetry.parity = 1; }
    /// Destructor
    virtual ~UniformField() {}
    
    vec3 B; ///< the uniform field vector
    
    /// Magnetic field at a specified point
    virtual vec3 fieldAt(const vec3&) const { return B; }
    /// Magnetic field averaged over specified line
    virtual vec3 fieldOverLine(Line) const { return B; }
    /// Magnetic field averaged over specified plane
    virtual vec3 fieldOverPlane(Plane) const { return B; }

    /// visualization
    virtual void _visualize() const {}
    
    /// Print info to stdout
    virtual void display() const { printf("[Uniform Field]\n"); }
};

#endif
