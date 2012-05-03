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
vec3 fieldLineIntegratorFunction(mdouble x, void* params) {
	fieldIntegratorParams* p = (fieldIntegratorParams*)params;
	vec3 pos = (p->l)->position(x);
	return (p->fs)->fieldAt(pos);
}

/// For integrating magnetic fields over a Plane using an Integrator
vec3 fieldPlaneIntegratorFunction(mdouble x, void* params) {
	fieldIntegratorParams p = *(fieldIntegratorParams*)params;
	Line l(p.p->position(x,-1.0),p.p->position(x,1.0));
	return (p.fs)->fieldOverLine(l);
}

//-------------------------------------

vec3 FieldSource::fieldOverLine(Line l) const {
	fieldIntegratorParams p;
	p.l = &l;
	p.fs = this;
	return lineIntegrator.integrate(&fieldLineIntegratorFunction,0.0,1.0,&p);
}

vec3 FieldSource::fieldOverPlane(Plane pl) const {
	fieldIntegratorParams p;
	p.p = &pl;
	p.fs = this;
	return planeIntegrator.integrate(&fieldPlaneIntegratorFunction,-1.0,1.0,&p)*0.5;
}
