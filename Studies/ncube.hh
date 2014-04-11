/* 
 * ncube.hh, part of the RotationShield program
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

#ifndef NCUBE_HH
/// Make sure this header is only loaded once
#define NCUBE_HH

#include "Visr.hh"
#include "Vec.hh"
#include "Vec_Null.hh"

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include <cassert>

/// N-dimensional rotation matrix
template<int N>
class NRotator {
public:
	/// constructor
	NRotator();
	/// destructor
	~NRotator() { gsl_matrix_free(M); }
	
	/// rotate by theta in the i-j plane
	void rotate(int i, int j, double th);
	/// apply rotation matrix to the given N-vector
	Vec<N,double> apply(const Vec<N,double>& v) const;
	
private:
	gsl_matrix* M; ///< the rotation matrix
};

template<int N>
NRotator<N>::NRotator() {
	M = gsl_matrix_calloc(N,N);
	for(int i=0; i<N; i++) gsl_matrix_set(M,i,i,1.0);
}

template<int N>
void NRotator<N>::rotate(int i, int j, double th) {
	gsl_matrix* R = gsl_matrix_calloc(N,N);
	for(int n=0; n<N; n++) gsl_matrix_set(R,n,n,1.0);
	double c = cos(th);
	double s = sin(th);
	gsl_matrix_set(R,i,i,c);
	gsl_matrix_set(R,i,j,-s);
	gsl_matrix_set(R,j,i,s);
	gsl_matrix_set(R,j,j,c);
	
	gsl_matrix* RM = gsl_matrix_alloc(N,N);
	assert(!gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1., R, M, 0., RM));
	gsl_matrix_free(R);
	gsl_matrix_free(M);
	M = RM;
}

template<int N>
Vec<N,double> NRotator<N>::apply(const Vec<N,double>& v) const{
	Vec<N,double> r = Vec<N,double>();
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			r[i] += gsl_matrix_get(M,i,j)*v[j];
	return r;
}



template<int N>
struct reprojectparams {
	void* p;
	Vec<3,double> (*f)(Vec<N,double>,void*);
	double d;
	int fixedplane;
	NRotator<N>* rotator;
};

template<int N>
vec3 reproject(Vec<N-1,double> v, void* params) {
	reprojectparams<N> p = *(reprojectparams<N>*)params;
	
	Vec<N,double> vN;
	for(int i=0; i<p.fixedplane; i++)
		vN[i] = v[i];
	vN[p.fixedplane] = p.d;
	for(int i=p.fixedplane+1; i<N; i++)
		vN[i] = v[i-1];
	return (*(p.f))(vN,p.p);
}

template<int N>
vec3 top_reproject(Vec<N-1,double> v, void* params) {
	reprojectparams<N> p = *(reprojectparams<N>*)params;
	
	Vec<N,double> vN;
	for(int i=0; i<p.fixedplane; i++)
		vN[i] = v[i];
	vN[p.fixedplane] = p.d;
	for(int i=p.fixedplane+1; i<N; i++)
		vN[i] = v[i-1];
	
	if(p.rotator) vN = p.rotator->apply(vN);
	
	vec3 v3;
	for(int i=0; i<3; i++) v3[i] = vN[i];
	return v3;
} 

/// class for drawing an N-dimensional hypercube
template<int N>
class NCube {
public:
	/// constructor
	NCube() {}
	/// destructor
	~NCube() {}
	/// draw the NCube
	void visualize(NRotator<N>* rotator = NULL, bool istop = true);
	/// draw the NCube as a face of a higher-dimensional NCube
	void sub_visualize(Vec<3,double> (*f)(Vec<N,double>,void*), void* params);
};


template<int N>
void NCube<N>::visualize(NRotator<N>* rotator, bool istop) {
	if(istop) {
		vsr::startRecording(true);
		vsr::clearWindow();
	}
	
	NCube<N-1> C;
	reprojectparams<N> p;
	p.rotator = rotator;
	for(int i=0; i<N; i++)
	{
		p.fixedplane = i;
		p.d = -0.5;
		C.sub_visualize(&top_reproject<N>,&p);
		p.d = 0.5;
		C.sub_visualize(&top_reproject<N>,&p);
	}
	
	if(istop) vsr::stopRecording();		
}

template<int N>
void NCube<N>::sub_visualize(Vec<3,double> (*f)(Vec<N,double>,void*), void* params) {
	NCube<N-1> C;
	
	reprojectparams<N> p;
	p.f = f;
	p.p = params;
	for(int fp=0; fp<N; fp++) {
		p.fixedplane = fp;
		p.d = -0.5;
		C.sub_visualize(&reproject<N>,(void*)&p);
		p.d = 0.5;
		C.sub_visualize(&reproject<N>,(void*)&p);
	}
}




template<>
void NCube<0>::sub_visualize(Vec<3,double> (*f)(Vec<0,double>,void*), void* params) {
	Vec<0,double> v0;
	Vec<3,double> v3 = (*f)(v0,params);
	vsr::setColor(1.0,0.0,0.0,1.0);
	vsr::dot(v3);
}

template<>
void NCube<1>::sub_visualize(Vec<3,double> (*f)(Vec<1,double>,void*), void* params) {
	NCube<0> C;
	
	reprojectparams<1> p;
	p.f = f;
	p.p = params;
	p.fixedplane = 0;
	p.d = -0.5;
	C.sub_visualize(&reproject<1>,(void*)&p);
	p.d = 0.5;
	C.sub_visualize(&reproject<1>,(void*)&p);
	
	Vec<1,double> v1;
	Vec<3,double> vs,ve;
	
	v1[0] = -0.5;
	vs = (*f)(v1,params);
	v1[0] = 0.5;
	ve = (*f)(v1,params);
	
	vsr::setColor(0.0,0.5,1.0,1.0);
	vsr::line(vs,ve);
}

template<>
void NCube<2>::sub_visualize(Vec<3,double> (*f)(Vec<2,double>,void*), void* params) {
	NCube<1> C;
	
	reprojectparams<2> p;
	p.f = f;
	p.p = params;
	for(int i=0; i<2; i++) {
		p.fixedplane = i;
		p.d = -0.5;
		C.sub_visualize(&reproject<2>,(void*)&p);
		p.d = 0.5;
		C.sub_visualize(&reproject<2>,(void*)&p);
	}
	
	Vec<2,double> v2;
	Vec<3,double> vc;
	float xyz[12];
	
	for(int dx = 0; dx<2; dx++) {
		for(int dy = 0; dy<2; dy++) {
			v2[0] = 0.98*(dx-0.5);
			v2[1] = 0.98*(dy-0.5);
			vc = (*f)(v2,params);
			xyz[3*(dx+2*dy) + 0] = vc[0];
			xyz[3*(dx+2*dy) + 1] = vc[1];
			xyz[3*(dx+2*dy) + 2] = vc[2];
		}
	}
	vsr::setColor(1.0,0.9,0.9,0.05);
	//glDisable(GL_DEPTH_TEST);
	vsr::filledquad(xyz);
	//glEnable(GL_DEPTH_TEST);
}


void cubeDemo()
{
	NCube<4> C;
	NRotator<4> R;
	while(1)
	{
		R.rotate(1,3,0.04);
		R.rotate(2,3,0.0314159265);
		R.rotate(0,3,.033737);
		//R.rotate(1,2,0.005);
		C.visualize(&R,true);
		usleep(20000);
	}
}

#endif
