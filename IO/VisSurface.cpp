/* 
 * VisSurface.cpp, part of the RotationShield program
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

#include "VisSurface.hh"
#include <cassert>

void VisSurface::_visualize() const {
	
	if(!mySurface) return;
	
	double dl1 = 1./vis_nx;
	double dl2 = 1./vis_ny;
	
	vsr::setColor(0.5,0.5,0.5,1);
	
	for(unsigned int nx = 0; nx < vis_nx; nx++) {
		for(unsigned int ny = 0; ny < vis_ny; ny++) {
		
			vec2 l0((nx+0.5)/vis_nx, (ny+0.5)/vis_ny);
			
			if(f_color) {
				vec4 clr = (*f_color)(l0,cparams);
				vsr::setColor(clr[0], clr[1] , clr[2], clr[3]);
			}
			
			float xyz[12];
			for(int dn1 = 0; dn1 < 2; dn1++) {
				for(int dn2 = 0; dn2 < 2; dn2++) {
					vec2 c = vec2(l0[0] + dl1*(dn1?0.5:-0.5), l0[1] + dl2*(dn2?0.5:-0.5));
					vec3 cx = (*mySurface)(c);
					if(f_height) {
						double h = (*f_height)(c,hparams);
						if(h)
							cx += mySurface->snorm(c,true)*h;
					}
					
					for(unsigned int i=0; i<3; i++)
						xyz[3*(2*dn1+dn2)+i] = cx[i];
				}
			}
			vsr::filledquad(xyz);
		}
	}
}

//--------------------------------------------

void SamplerSurface::clear_data() {
	for(unsigned int i=0; i<G.size(); i++)
		delete G[i];
	G.clear();
}

void SamplerSurface::make_grids(unsigned int nx, unsigned int ny, unsigned int ng) {
	clear_data();
	NX = nx;
	NY = ny;
	for(unsigned int i=0; i<ng; i++)
		G.push_back(new BicubicGrid(NX,NY));
}

void SamplerSurface::sample_func(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* fparams) {
	assert(f);
	for(unsigned int x=0; x<NX; x++) {
		for(unsigned int y=0; y<NY; y++) {
			vec2 l(float(x)/(NX-1),float(y)/(NY-1));
			vec2 v(ll[0]*(1-l[0])+ur[0]*l[0], ll[1]*(1-l[1])+ur[1]*l[1]);
			mvec z = (*f)(ll,fparams);
			for(unsigned int i=0; i<z.size(); i++) {
				if(i>G.size()) break;
				G[i]->set(x,y,z[i]);
			}
		}
	}
}

mvec SamplerSurface::eval(vec2 l) const {
	mvec v(G.size());
	for(unsigned int i=0; i<G.size(); i++)
		v[i] = (*G[i])(l[0],l[1]);
	return v;
}

//--------------------------------------------

double h_fromSampler1(vec2 l, const void* smplr) {
	const SamplerSurface* S = (const SamplerSurface*)smplr;
	return S->eval(0,l);
}

double h_fromSamplerMag(vec2 l, const void* smplr) {
	const SamplerSurface* S = (const SamplerSurface*)smplr;
	return S->eval(l).mag();
}

vec4 c_fromSampler2(vec2 l, const void* smplr) {
	const SamplerSurface* S = (const SamplerSurface*)smplr;
	mvec v = S->eval(l);
	vec3 hsv = vec3( atan2(v[1],v[0]), 1.0, 1.0 );
	vec3 rgb = hsv2rgb(hsv);
	return vec4(rgb[0],rgb[1],rgb[2],1);
}

//---------------------------------------------

void vis_integrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* params) {

	SamplerSurface SS;
	SS.make_grids(20,20,3);
	SS.sample_func(f,ll,ur,params);
	
	for(unsigned int i=0; i<SS.G.size(); i++) {
		SS.G[i]->scale_zrange(-0.2,0.2);
	}
	
	VisSurface VS(new Plane3D);
	VS.hparams = VS.cparams = &SS;
	VS.f_height = &h_fromSamplerMag;
	VS.f_color = &c_fromSampler2;
	
	VS.visualize();
}
