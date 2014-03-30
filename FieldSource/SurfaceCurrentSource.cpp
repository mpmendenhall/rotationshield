/* 
 * SurfaceCurrentSource.cpp, part of the RotationShield program
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

#include "SurfaceCurrentSource.hh"
#include "Color.hh"
#include <cmath>
#include <algorithm>

vec3 SurfaceCurrentSource::dI_contrib(const vec2& l) const {
	assert(mySurface && sj);
	
	// surface current element dl
	vec2 sdl = (*sj)(l, sjparams);
	vec3 dl0 = mySurface->deriv(l,0);
	vec3 dl1 = mySurface->deriv(l,1);
	double ml0 = dl0.mag();
	double ml1 = dl1.mag();
	if(!ml0 || !ml1) return vec3(0,0,0);
	
	double dA = cross(dl0,dl1).mag();
	return vec3(dl0/ml0 * sdl[0] + dl1/ml1 * sdl[1])*dA;
}

vec3 SurfaceCurrentSource::fieldAt_contrib_from(const vec3& v, const vec2& l) const {
	assert(mySurface);
	
	vec3 x0 = (*mySurface)(l);
	vec3 dI = dI_contrib(l);
		
	// Biot-Savart law
	vec3 r = v - x0;
	double mr = r.mag();
	if(!mr) return vec3(0,0,0);
	vec3 B = cross(dI,r)/(4.*M_PI*mr*mr*mr);
	if(!(B[0]==B[0] && B[1]==B[1] && B[2]==B[2])) {
		std::cout << v << l << mr << B << std::endl;
		assert(false);
	}
	return B;
}

void SurfaceCurrentSource::_visualize() const {
	
	if(!mySurface || !sj) return;
	
	double dl1 = 1./vis_n1;
	double dl2 = 1./vis_n2;
	double logmax = 3.0;
	
	for(unsigned int n1 = 0; n1 < vis_n1; n1++) {
		for(unsigned int n2 = 0; n2 < vis_n2; n2++) {
		
			vec2 l0((n1+0.5)/vis_n1, (n2+0.5)/vis_n2);
			
			vec2 j = (*sj)(l0,sjparams);
			
			float smag = j.mag2();
			if(smag)
				smag = std::max(0.,std::min(1.0,0.1*(log(smag)+10-logmax)));
			vec3 hsv = vec3( atan2(j[1],j[0]), smag, 1.0 );
			vec3 rgb = hsv2rgb(hsv);
			vsr::setColor(rgb[0], rgb[1] , rgb[2], 1.0);
		
			float xyz[12];
			for(int dn1 = 0; dn1 < 2; dn1++) {
				for(int dn2 = 0; dn2 < 2; dn2++) {
					vec2 c = vec2(l0[0] + dl1*(dn1?0.5:-0.5), l0[1] + dl2*(dn2?0.5:-0.5));
					vec3 cx = (*mySurface)(c);
					for(unsigned int i=0; i<3; i++)
						xyz[3*(2*dn1+dn2) + i] = cx[i];
				}
			}
			vsr::filledquad(xyz);
		}
	}
}

void SurfaceCurrentSource::vis_coords(const vec2& l, double s) const {
	if(!mySurface) return;
	
	vec3 o = (*mySurface)(l);
	vec3 dx = mySurface->deriv(l,0).normalized();
	vec3 dy = mySurface->deriv(l,1).normalized();
	vec3 dz = cross(dx,dy);
	
	vsr::line(o, o+dx*s);
	vsr::line(o, o+dy*s);
	vsr::line(o, o+dz*s);
}

mvec J_dA(vec2 l, void* params) {
	SurfaceCurrentSource& S = *(SurfaceCurrentSource*)params;
	return mvec(S.dI_contrib(l));
}

vec3 SurfaceCurrentSource::netCurrent(vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	assert(mySurface);
	myIntegrator.setMethod(INTEG_GSL_QAG);
	mvec J = subdividedIntegral(&J_dA, 3, (void*)this, ll, ur, ndx, ndy);
	return vec3(J[0],J[1],J[2]);
}



