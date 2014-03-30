/* 
 * FieldSource.cpp, part of the RotationShield program
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

#include "FieldSource.hh"

Integrator FieldSource::lineIntegrator = Integrator();
Integrator FieldSource::planeIntegrator = Integrator();

/// Used for numerical integration of magnetic fields over Lines/Planes
struct fieldIntegratorParams {
	const Line* l; //< Line over which to integrate
	const Plane* p; //< Plane over which to integrate
	const FieldSource* fs; //< source of the fields being integrated
};

/// For integrating magnetic fields over a Line using an Integrator
mvec fieldLineIntegratorFunction(double x, void* params) {
	fieldIntegratorParams* p = (fieldIntegratorParams*)params;
	vec3 pos = (p->l)->position(x);
	return mvec((p->fs)->fieldAt(pos));
}

/// For integrating magnetic fields over a Plane using an Integrator
mvec fieldPlaneIntegratorFunction(double x, void* params) {
	fieldIntegratorParams p = *(fieldIntegratorParams*)params;
	Line l(p.p->position(x,-1.0),p.p->position(x,1.0));
	return mvec((p.fs)->fieldOverLine(l));
}

//-------------------------------------

vec3 FieldSource::fieldOverLine(Line l) const {
	fieldIntegratorParams p;
	p.l = &l;
	p.fs = this;
	mvec v = lineIntegrator.integrate(&fieldLineIntegratorFunction,0.0,1.0,&p);
	return vec3(v[0],v[1],v[2]);
}

vec3 FieldSource::fieldOverPlane(Plane pl) const {
	fieldIntegratorParams p;
	p.p = &pl;
	p.fs = this;
	mvec v = planeIntegrator.integrate(&fieldPlaneIntegratorFunction,-1.0,1.0,&p)*0.5;
	return vec3(v[0],v[1],v[2]);
}
