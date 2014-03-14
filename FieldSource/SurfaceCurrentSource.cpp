#include "SurfaceCurrentSource.hh"
#include "Color.hh"
#include <cmath>

vec3 SurfaceCurrentSource::dI_contrib(const vec2& l, vec3& xout) const {
	assert(mySurface && sj);
	
	// position of surface element
	xout = (*mySurface)(l);

	// surface current element dl
	vec2 sdl = (*sj)(l, sjparams);
	vec3 dl0 = mySurface->deriv(l,0);
	vec3 dl1 = mySurface->deriv(l,1);
	double ml0 = dl0.mag();
	double ml1 = dl1.mag();
	double dA = cross(dl0,dl1).mag();
	
	return vec3(dl0/ml0 * sdl[0] + dl1/ml1 * sdl[1])*dA;
}

vec3 SurfaceCurrentSource::fieldAt_contrib_from(const vec3& v, const vec2& l) const {
	vec3 x0;
	vec3 dl = dI_contrib(l,x0);
		
	// Biot-Savart law, transitioning to near-field constant field approximation TODO
	vec3 r = v - x0;
	double mr = r.mag();
	
	if(mr < 1e-9) return vec3(0,0,0);	//< zero out singularity
	
	return cross(dl,r)/(4.*M_PI*mr*mr*mr);
	
	
	//if(mr > 0.05)
	//if(mr > 0.02) return vec3();
	//return vec3();
	
}

void SurfaceCurrentSource::_visualize() const {
	
	if(!mySurface || !sj) return;
	
	mdouble dl1 = 1./vis_n1;
	mdouble dl2 = 1./vis_n2;
	mdouble logmax = 3.0;
	
	for(unsigned int n1 = 0; n1 < vis_n1; n1++) {
		for(unsigned int n2 = 0; n2 < vis_n2; n2++) {
		
			vec2 l0((n1+0.5)/vis_n1, (n2+0.5)/vis_n2);
			
			vec2 j = (*sj)(l0,sjparams);
			
			vec3 hsv = vec3( atan2(j[1],j[0]), 1.0, 1.0 );
			vec3 rgb = hsv2rgb(hsv);
			float smag = j.mag2();
			if(smag)
				smag = max(0,min(1.0,0.1*(log(smag)+10-logmax)));
			vsr::setColor(rgb[0], rgb[1] , rgb[2], 1.0*smag);
		
			float xyz[12];
			for(int dn1 = 0; dn1 < 2; dn1++) {
				for(int dn2 = 0; dn2 < 2; dn2++) {
					vec2 c = vec2(l0[0] + dl1*(dn1?0.5:-0.5), l0[1] + dl2*(dn2?0.5:-0.5));
					vec3 cx = (*mySurface)(c);
					for(unsigned int i=0; i<3; i++)
						xyz[3*(2*dn1+dn2)+i] = cx[i];
				}
			}
			vsr::filledquad(xyz);
		}
	}
}

void SurfaceCurrentSource::vis_coords(const vec2& l, double s) const {
	vec3 o = (*mySurface)(l);
	vec3 dx = mySurface->deriv(l,0).normalized();
	vec3 dy = mySurface->deriv(l,1).normalized();
	vec3 dz = cross(dx,dy);
	
	vsr::line(o, o+dx*s);
	vsr::line(o, o+dy*s);
	vsr::line(o, o+dz*s);
}

