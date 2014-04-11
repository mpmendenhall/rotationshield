/* 
 * PlaneSource.hh, part of the RotationShield program
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

#ifndef PLANESOURCE_HH
/// Make sure this header is only loaded once
#define PLANESOURCE_HH

#include "PlanarElement.hh"
#include "Matrix.hh"

/// Magnetic field source due to a plane of current (e.g. a bound surface current on a magnetic shield)
class PlaneSource: public PlanarElement {
public:
	/// Constructor
	PlaneSource(Plane cp, double mu): PlanarElement(cp), murel(mu) { setRmat(); }
	
	/// number of degrees of freedom
	virtual unsigned int nDF() const { return 2; }
	
	/// state in response to applied field
	virtual mvec responseToFieldSource(const FieldSource* f) const { return rmat * f->fieldOverPlane(p); }
	
	/// interaction matrix with another ReactiveElement
	virtual mmat interactionWithComponents(const ReactiveElement* e) const { return rmat * e->fieldAtComponents(p.o); }
	/// interaction matrix with another ReactiveElement
	virtual mvec interactionWith(const ReactiveElement* e) const { return rmat * e->fieldAt(p.o); }

	/// replicate around a new plane
	virtual PlanarElement* replicate(Plane pl) const { return new PlaneSource(pl,murel); }
	
	/// generate initial reference element
	virtual PlanarElement* reference(annulusSpec a) const { return new PlaneSource(Plane(a),murel); }
	
	/// replicate reference element to other angle
	virtual PlanarElement* replicateRotated(double th) const { return new PlaneSource(p.zrotated(th),murel); }

	/// field components due to each DF at given point
	virtual mmat fieldAtComponents(vec3 p0) const;
			
private:
	double murel;	///< relative permeability
	mmat rmat;		///< response matrix to applied field
	
	/// generate correct response matrix to applied fields
	void setRmat() {
		Matrix<2,3,double> M = Matrix<2,3,double>();
		M(0,1) = M(1,0) = 2.0*(1.0-murel)/(1.0+murel);
		rmat = mmat(M*p.projectionMatrix());		
	}
	
};

#endif
