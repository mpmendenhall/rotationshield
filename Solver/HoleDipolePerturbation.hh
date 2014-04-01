/* 
 * HoleDipolePerturbation.hh, part of the RotationShield program
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

#ifndef HOLEDIPOLEPERTURBATION_HH
#define HOLEDIPOLEPERTURBATION_HH 1

#include "DipoleSource.hh"
#include "ReactiveSet.hh"
#include "SurfaceGeometry.hh"
#include "MagRS.hh"

/// A small hole in a superconducting sheet
class HoleDipolePerturbation: public ReactiveSet, public DipoleSource {
public:
	/// constructor
	HoleDipolePerturbation(const SurfaceGeometry& SG, vec2 l, double r): ReactiveSet(1), DipoleSource(SG(l),vec3()), a(r), dh(0.05), mySurface(SG), surfacePos(l) {}
	
	//===================================== ReactiveSet subclass
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return 2; }
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0);
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip);
	//=====================================
	
	/// Visualize the interactor
	virtual void _visualize() const;
	
	double a;	//< hole radius
	double dh;	//< height above surface to evaluate response field
	
protected:

	//===================================== ReactiveSet subclass
	/// called when a DF is set
	virtual void _setDF(unsigned int DF, double v);
	//=====================================
	
	const SurfaceGeometry& mySurface;	//< surface to which this is "attached"
	vec2 surfacePos;					//< position on surface
};

#endif
