#include "SurfaceSource.hh"
#include <iostream>

struct SurfaceSourceIntegParams {
	const SurfaceSource* S;
	vec3 v;
};

mvec SSdA(mdouble x, mdouble y, void* params) {
	SurfaceSourceIntegParams& p = *(SurfaceSourceIntegParams*)params;
	return mvec(p.S->fieldAt_contrib_from(p.v,x,y));
}

mdouble clamp(mdouble m) { return m<0? 0 : m<1? m : 1; }

vec3 SurfaceSource::fieldAt(const vec3& v, const vec2& ll, const vec2& ur) const {
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	mvec B = myIntegrator.integrate(&SSdA, clamp(ll[0]), clamp(ur[0]), clamp(ll[1]), clamp(ur[1]), &p);
	return vec3(B[0],B[1],B[2]);
}

void SurfaceSource::displayContribGrid(const vec3& v, unsigned int nx, unsigned int ny) const {
	std::cout << "Field contributions to " << v << std::endl;
	vec3 B;
	for(unsigned int ix=0; ix<nx; ix++) {
		mdouble x = float(ix)/(nx-1);
		std::cout << x << "\t";
		for(unsigned int iy=0; iy<ny; iy++) {
			mdouble y = float(iy)/(ny-1);
			vec3 Bi = fieldAt_contrib_from(v,x,y);
			std::cout << "\t" << Bi;
			B += Bi;
		}
		std::cout << std::endl;
	}
	B /= nx*ny;
	std::cout << "Total field: " << B << std::endl;
}
