/* 
 * SurfaceCurrentRS.hh, part of the RotationShield program
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

/// \file SurfaceCurrentRS.hh \brief ReactiveSet for a surface current responding to magnetic field conditions

#ifndef SURFACECURRENTRS_HH
/// Makes sure to only load this file once
#define SURFACECURRENTRS_HH 1

#include "SurfaceCurrentSource.hh"
#include "InterpolatingRS.hh"
#include "MagRS.hh"

/// defines surface boundary condition by producing surface current in response to applied magnetic field
class SurfaceI_Response {
public:
	/// constructor
	SurfaceI_Response(double mu = 0) { setMu(mu); }
	
	/// generate correct response matrix to applied fields for given relative permeability
	void setMu(double mu) {
		murel = mu;
		rmat2 = Matrix<2,3,double>();
		rmat2(0,1) = 2*(murel-1)/(murel+1);
		rmat2(1,0) = -rmat2(0,1);
	}
	
	double murel;					//< relative permeability
	Matrix<2,3,double> rmat2;		//< response matrix to applied field
};


/// Continuous surface current responding to magnetic field
class SurfaceCurrentRS: public MagF_Responder, public SurfaceCurrentSource, public InterpolatingRS2D {
public:
	/// constructor
	SurfaceCurrentRS(SurfaceGeometry* SG, unsigned int nph, unsigned int nz, const std::string& nm = "SurfaceCurrentRS");

	// ReactiveSet subclassed functions
	//=====================================
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0);
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip);
	//=====================================
	
	/// calculate response to incident field
	virtual void calculateIncident(const FieldSource& f);
	
	/// set surface response at all points
	void setSurfaceResponse(SurfaceI_Response r);
	
	/// visualization routine
	virtual void _visualize() const;
	/// visualize current vectors
	virtual void vis_i_vectors(double s = 0.5, double mx = 0.2) const;
	
	/// get center surface coordinates for i^th element
	vec2 surf_coords(unsigned int i) const { return vec2( ((i/nPhi)+0.5)/nZ, ((i%nPhi)+0.5)/nPhi); }
	/// get surface coordinate range for element (not including extended interpolation effects)
	void element_surface_range(unsigned int i, vec2& ll, vec2& ur) const;
	
	/// get interpolated surface current response
	vec2 eval_J(const vec2& p) const { return vec2( (*G[0])(p[0],p[1]), (*G[1])(p[0],p[1]) ); }
	
	/// set up current loop response around given ring of z elements
	void set_current_loop(unsigned int z, double i=1.0, bool phidir = true);
	
protected:
	
	/// one surface element's reaction
	mvec subelReaction(ReactiveSet* R);
	
	unsigned int ixn_el;	//< interacting element currently being probed
	
	std::vector<SurfaceI_Response> sdefs;	//< surface response definitions at each site
};

#endif
