#include "SurfaceSource.hh"
#include "VisSurface.hh"
#include <iostream>
#include <cassert>
#include <algorithm>

struct SurfaceSourceIntegParams {
	const SurfaceSource* S;
	vec3 v;
	const Matrix<2,3,double>* M2;
	const Matrix<3,3,double>* M3;
};

mvec SSdA(vec2 l, void* params) {
	SurfaceSourceIntegParams& p = *(SurfaceSourceIntegParams*)params;
	vec3 B = p.S->fieldAt_contrib_from(p.v,l);
	if(p.M2) return mvec( (*p.M2)*B );
	if(p.M3) return mvec( (*p.M3)*B );
	return mvec(B);
}

vec3 SurfaceSource::fieldAt(const vec3& v, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	p.M2 = NULL;
	p.M3 = NULL;
	
	mvec B = mySurface->subdividedIntegral(&SSdA, 3, &p, ll, ur, ndx, ndy);
	assert(B.size()==3);
	return vec3(B[0],B[1],B[2]);

}

vec2 SurfaceSource::fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	p.M2 = &M;
	p.M3 = NULL;
	
	mvec MB = mySurface->subdividedIntegral(&SSdA, 2, &p, ll, ur, ndx, ndy);
	assert(MB.size()==2);
	return vec2(MB[0],MB[1]);
}

vec3 SurfaceSource::fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	
	SurfaceSourceIntegParams p;
	p.S = this;
	p.v = v;
	p.M2 = NULL;
	p.M3 = &M;
	
	mvec MB = mySurface->subdividedIntegral(&SSdA, 3, &p, ll, ur, ndx, ndy);
	assert(MB.size()==3);
	return vec3(MB[0],MB[1],MB[2]);
}

//---------


void SurfaceSource::displayContribGrid(const vec3& v, unsigned int nx, unsigned int ny) const {
	std::cout << "Field contributions to " << v << std::endl;
	vec3 B;
	for(unsigned int ix=0; ix<nx; ix++) {
		double x = float(ix)/(nx-1);
		std::cout << x << "\t";
		for(unsigned int iy=0; iy<ny; iy++) {
			double y = float(iy)/(ny-1);
			vec3 Bi = fieldAt_contrib_from(v,vec2(x,y));
			std::cout << "\t" << Bi;
			B += Bi;
		}
		std::cout << std::endl;
	}
	B /= nx*ny;
	std::cout << "Total field: " << B << std::endl;
}

void SurfaceSource::visualize_line(vec2 s, vec2 e) const {
	if(!vis_on) return;
	unsigned int nseg = std::max((unsigned int)std::max(fabs(s[0]-e[0])*vis_n1, fabs(s[1]-e[1])*vis_n2),(unsigned int)1);
	for(unsigned int i=0; i<nseg; i++) {
		double l1 = double(i)/nseg;
		double l2 = double(i+1)/nseg;
		vsr::line((*mySurface)(s*(1-l1)+e*l1), (*mySurface)(s*(1-l2)+e*l2));
	}
}

void SurfaceSource::visualize_region(vec2 ll, vec2 ur) const {
	visualize_line(ll,vec2(ur[0],ll[1]));
	visualize_line(vec2(ur[0],ll[1]),ur);
	visualize_line(ur,vec2(ll[0],ur[1]));
	visualize_line(vec2(ll[0],ur[1]),ll);
}



