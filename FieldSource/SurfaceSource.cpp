#include "SurfaceSource.hh"
#include <iostream>

struct SurfaceSourceIntegParams {
	const SurfaceSource* S;
	vec3 v;
};

mvec SSdA(vec2 l, void* params) {
	SurfaceSourceIntegParams& p = *(SurfaceSourceIntegParams*)params;
	return mvec(p.S->fieldAt_contrib_from(p.v,l));
}

mdouble clamp(mdouble m) { return m<0? 0 : m<1? m : 1; }

vec3 SurfaceSource::fieldAt(const vec3& v, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	
	ndx = ndx?ndx:dflt_integrator_ndivs_x;
	ndy = ndy?ndx:dflt_integrator_ndivs_y;
	for(unsigned int i=0; i<2; i++) {
		ll[i] = clamp(ll[i]);
		ur[i] = (clamp(ur[i])-ll[i])/(i?ndy:ndx);
	}
	
	vec3 B;
	for(unsigned int nx=0; nx<ndx; nx++) {
		for(unsigned int ny=0; ny<ndy; ny++) {
			vec2 lll =ll + vec2(nx*ur[0],ny*ur[1]);
			mvec Bi = myIntegrator.integrate2D(&SSdA, lll, lll+ur, &p);
			B += vec3(Bi[0],Bi[1],Bi[2]);
		}
	}
	return B;
}

void SurfaceSource::displayContribGrid(const vec3& v, unsigned int nx, unsigned int ny) const {
	std::cout << "Field contributions to " << v << std::endl;
	vec3 B;
	for(unsigned int ix=0; ix<nx; ix++) {
		mdouble x = float(ix)/(nx-1);
		std::cout << x << "\t";
		for(unsigned int iy=0; iy<ny; iy++) {
			mdouble y = float(iy)/(ny-1);
			vec3 Bi = fieldAt_contrib_from(v,vec2(x,y));
			std::cout << "\t" << Bi;
			B += Bi;
		}
		std::cout << std::endl;
	}
	B /= nx*ny;
	std::cout << "Total field: " << B << std::endl;
}
