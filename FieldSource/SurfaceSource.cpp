#include "SurfaceSource.hh"

struct SurfaceSourceIntegParams {
	const SurfaceSource* S;
	vec3 v;
};

vec3 SSdA(mdouble x, mdouble y, void* params) {
	SurfaceSourceIntegParams& p = *(SurfaceSourceIntegParams*)params;
	return p.S->fieldAt_contrib_from(p.v,x,y);
}

vec3 SurfaceSource::fieldAt(const vec3& v) const {
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	return myIntegrator.integrate(&SSdA, 0, 1, 0, 1, &p);
	
}
