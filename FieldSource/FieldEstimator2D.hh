/* 
 * FieldEstimator2D.hh, part of the RotationShield program
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

/// \file FieldEstimator2D.hh \brief 2D magnetic fields

#ifndef FIELDESTIMATOR2D_HH
/// Make sure this header is only loaded once
#define FIELDESTIMATOR2D_HH

#include "FieldSource.hh"
#include <vector>

/// Rough estimate of field magnitude (in 2D plane) due to wires perpendicular to shield, used to optimize shield gridding for cos theta coil endcaps
class FieldEstimator2D {
public:
    /// Constructor
    FieldEstimator2D() {}
    /// Destructor
    virtual ~FieldEstimator2D() {}
    /// Estimated field at a point
    virtual vec2 estimateAt(const vec2& v) const;
    /// estimate rate of change in given direction
    virtual vec2 derivAt(const vec2& v, const vec2& dx) const;
    /// Add a "line source" perpendicular to plane of interest
    void addsource(const vec2& v, double j) { sources.push_back(v); currents.push_back(j); }
protected:
    std::vector<vec2> sources;
    std::vector<double> currents;
};

/// Field estimator based off 3D field source
class FieldEstimator2Dfrom3D: public FieldEstimator2D {
public:
    /// Constructor
    FieldEstimator2Dfrom3D(FieldSource* FS): myFS(FS) { assert(myFS); myFS->retain(); }
    /// Destructor
    virtual ~FieldEstimator2Dfrom3D() { myFS->release(); }
    /// Estimated field at a point
    virtual vec2 estimateAt(const vec2& v) const;
protected:
    FieldSource* myFS;
};

#endif
