/* 
 * MagRS.hh, part of the RotationShield program
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

#ifndef MAGRS_HH
#define MAGRS_HH 1

#include "ReactiveSet.hh"
#include "FieldSource.hh"
#include "MixedSource.hh"

/// Collections of magnetic-field-responding ReactiveSets
class MagRSCombiner: public ReactiveSetCombiner, public FieldSource {
public:
	/// constructor
	MagRSCombiner(unsigned int nphi): ReactiveSetCombiner(nphi) {}
	
	/// add a ReactiveSet
	virtual void addSet(ReactiveSet* R);
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const;
	/// visualization reoutine
	virtual void _visualize() const;
};

/// Applied magnetic field source
class MagExtField: public ReactiveSet, public MixedSource {
public:
	/// constructor
	MagExtField(): ReactiveSet(0), MixedSource() {}

	//=====================================
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return 0; }
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet*, unsigned int = 0) { return mvec(); }
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip);
	//=====================================
	
protected:

	//=====================================
	/// called when a DF is set
	virtual void _setDF(unsigned int, double) { }
	/// optional routine for setting entire state vector at once
	virtual void _setDFv(const mvec&) { }
	//=====================================
};

/// Magnetic field interaction protocol class singleton
class BField_Protocol {
public:
	vec3 x;	//< position
	vec3 B;	//< magnetic field
	const Matrix<2,3,double>* M2;	//< optional 3-to-2 transform matrix
	const Matrix<3,3,double>* M3;	//< optional 3-to-3 transform matrix
	vec2 M2B;						//< transformed field
	const void* caller;				//< pointer to caller
	static BField_Protocol* BFP;	//< instance to use
};

#endif
