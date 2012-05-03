#ifndef NCUBE_HH
#define NCUBE_HH 1

#include "Visr.hh"
#include "Vec.hh"
#include "Vec_Null.hh"

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

/// N-dimensional rotation matrix
template<int N>
class NRotator {
public:
	/// constructor
	NRotator(): M(N,N) { for(int i=0; i<N; i++) M[i][i] = 1.0; }
	/// destructor
	~NRotator() {}
	
	/// rotate by theta in the i-j plane
	void rotate(int i, int j, mdouble th);
	/// apply rotation matrix to the given N-vector
	Vec<N,mdouble> apply(const Vec<N,mdouble>& v) const;
	
private:
	CLHEP::HepMatrix M; //< the rotation matrix
};

template<int N>
void NRotator<N>::rotate(int i, int j, mdouble th) {
	CLHEP::HepMatrix R(N,N);
	for(int n=0; n<N; n++) R[n][n] = 1.0;
	mdouble c = cos(th);
	mdouble s = sin(th);
	R[i][i] = c;
	R[i][j] = -s;
	R[j][i] = s;
	R[j][j] = c;
	
	M = R*M;
}

template<int N>
Vec<N,mdouble> NRotator<N>::apply(const Vec<N,mdouble>& v) const{
	Vec<N,mdouble> r = Vec<N,mdouble>();
	for(int i=0; i<N; i++)
		for(int j=0; j<N; j++)
			r[i] += M[i][j]*v[j];
	return r;
}



template<int N>
struct reprojectparams {
	void* p;
	Vec<3,mdouble> (*f)(Vec<N,mdouble>,void*);
	mdouble d;
	int fixedplane;
	NRotator<N>* rotator;
};

template<int N>
vec3 reproject(Vec<N-1,mdouble> v, void* params) {
	reprojectparams<N> p = *(reprojectparams<N>*)params;
	
	Vec<N,mdouble> vN;
	for(int i=0; i<p.fixedplane; i++)
		vN[i] = v[i];
	vN[p.fixedplane] = p.d;
	for(int i=p.fixedplane+1; i<N; i++)
		vN[i] = v[i-1];
	return (*(p.f))(vN,p.p);
}

template<int N>
vec3 top_reproject(Vec<N-1,mdouble> v, void* params) {
	reprojectparams<N> p = *(reprojectparams<N>*)params;
	
	Vec<N,mdouble> vN;
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
	void sub_visualize(Vec<3,mdouble> (*f)(Vec<N,mdouble>,void*), void* params);
};


template<int N>
void NCube<N>::visualize(NRotator<N>* rotator, bool istop) {
	if(istop) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); }
	
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
	
	if(istop) vsr::Visr::W->stopRecording();		
}

template<int N>
void NCube<N>::sub_visualize(Vec<3,mdouble> (*f)(Vec<N,mdouble>,void*), void* params) {
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
void NCube<0>::sub_visualize(Vec<3,mdouble> (*f)(Vec<0,mdouble>,void*), void* params) {
	Vec<0,mdouble> v0;
	Vec<3,mdouble> v3 = (*f)(v0,params);
	vsr::Visr::W->setColor(1.0,0.0,0.0,1.0);
	vsr::Visr::W->dot(v3);
}

template<>
void NCube<1>::sub_visualize(Vec<3,mdouble> (*f)(Vec<1,mdouble>,void*), void* params) {
	NCube<0> C;
	
	reprojectparams<1> p;
	p.f = f;
	p.p = params;
	p.fixedplane = 0;
	p.d = -0.5;
	C.sub_visualize(&reproject<1>,(void*)&p);
	p.d = 0.5;
	C.sub_visualize(&reproject<1>,(void*)&p);
	
	Vec<1,mdouble> v1;
	Vec<3,mdouble> vs,ve;
	
	v1[0] = -0.5;
	vs = (*f)(v1,params);
	v1[0] = 0.5;
	ve = (*f)(v1,params);
	
	vsr::Visr::W->setColor(0.0,0.5,1.0,1.0);
	vsr::Visr::W->line(vs,ve);
}

template<>
void NCube<2>::sub_visualize(Vec<3,mdouble> (*f)(Vec<2,mdouble>,void*), void* params) {
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
	
	Vec<2,mdouble> v2;
	Vec<3,mdouble> vc;
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
	vsr::Visr::W->setColor(1.0,0.9,0.9,0.05);
	//glDisable(GL_DEPTH_TEST);
	vsr::Visr::W->filledquad(xyz);
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
		usleep(50000);
	}
}

#endif
