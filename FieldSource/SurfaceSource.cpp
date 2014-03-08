#include "SurfaceSource.hh"

struct SurfaceSourceIntegParams {
	const SurfaceSource* S;
	vec3 v;
};

mvec SSdA(mdouble x, mdouble y, void* params) {
	SurfaceSourceIntegParams& p = *(SurfaceSourceIntegParams*)params;
	return mvec(p.S->fieldAt_contrib_from(p.v,x,y));
}

vec3 SurfaceSource::fieldAt(const vec3& v) const {
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	mvec B = myIntegrator.integrate(&SSdA, 0, 1, 0, 1, &p);
	return vec3(B[0],B[1],B[2]);
	
}
