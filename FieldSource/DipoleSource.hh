/* 
 * DipoleSource.hh, part of the RotationShield program
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

#ifndef DIPOLESOURCE_HH
/// Make sure this header is only loaded once
#define DIPOLESOURCE_HH

#include "FieldSource.hh"


/// Magenetic dipole field source
class DipoleSource: public FieldSource {
public:
    /// Constructor
    DipoleSource(vec3 X, vec3 M): FieldSource("DipoleSource"), x(X), m(M) { }
    
    /// Field at a specified point
    virtual vec3 fieldAt(const vec3& v) const;
    /// Visualize the field source
    virtual void _visualize() const;
    
protected:
    vec3 x;     ///< location
    vec3 m;     ///< dipole moment
};

#endif
