/* 
 * SurfaceCurrentSource.hh, part of the RotationShield program
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

/// \file SurfaceCurrentSource.hh \brief Magnetic field from a surface current density

#ifndef SURFACECURRENTSOURCE_HH
/// Makes sure to only load this file once
#define SURFACECURRENTSOURCE_HH 1

#include "SurfaceSource.hh"

class SurfaceCurrentSource: public SurfaceSource {
public:
	/// constructor
	SurfaceCurrentSource(SurfaceGeometry* SG = NULL, const std::string& nm = "SurfaceCurrentSource"):
		SurfaceSource(SG,nm), sj(NULL), sjparams(NULL) {}
		
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, const vec2& l) const;

	/// net current from specified region
	virtual vec3 netCurrent(vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	/// Net current of the FieldSource
	virtual vec3 net_current() const { return netCurrent(vec2(0,0),vec2(1,1)); }
	
	vec2 (*sj)(vec2, void*);	//< surface current density as a function of surface coordinate
	void* sjparams;				//< extra parameters for surface function
	
	/// Visualize the field source
	virtual void _visualize() const;
	
	/// visualize local coordinate axes
	virtual void vis_coords(const vec2& l, double s = 0.02) const;

	/// current element dI from surface coordinate l
	vec3 dI_contrib(const vec2& l) const;
	
};


#endif
