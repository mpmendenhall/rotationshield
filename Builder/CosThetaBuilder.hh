#ifndef COSTHETABUILDER_HH
#define COSTHETABUILDER_HH 1


#include "MixedSource.hh"
#include <vector>

/// Cos theta coil with fourier-series of wire shifts from normal
mdouble shiftPositioner(unsigned int i, unsigned int ncoils, void* params) {
	i = i%ncoils;
	mdouble x0 = (0.5+i)/(mdouble)ncoils;
	mdouble x = x0;
	VarVec<mdouble>* shift = (VarVec<mdouble>*)params;
	if(shift)
		for(unsigned int k=0; k<shift->size(); k++) x += (*shift)[k]*sin(PI*(mdouble)(k+1)*x0);
	mdouble y = sqrt(1-x*x);
	return atan2(y,x);
}

/// Cos theta coil positioner using "k" parameter
mdouble alarconKPositioner(unsigned int i, unsigned int ncoils, void* params) {	
	i = i%ncoils;
	mdouble k = *(mdouble*)params;
	mdouble x0 = (0.5+i)/(mdouble)ncoils;
	return atan2(sqrt(1/(x0*x0)-1),(1-k)*(1-k));
}

vec3 simpletrans(int N, int i, int a1, int a2, int a3, void* params) {
	if(params)
		return *(vec3*)params;
	return vec3();
}

struct doubleTrans {
	vec3 trans1;
	vec3 trans2;
};

/// different translation for +z,-z ends
vec3 fancytrans(int N, int i, int a1, int a2, int a3, void* params) {
	if(!params)
		return vec3();
	doubleTrans* t = (doubleTrans*)params;
	if(a1>0)
		return t->trans1;
	return t->trans2;
}


class CosThetaBuilder {
public:
	CosThetaBuilder(unsigned int n, mdouble r, mdouble l,
					mdouble (*afunc)(unsigned int, unsigned int, void*) = &shiftPositioner,
					vec3 (*tfunc)(int, int, int, int, int, void*) = &simpletrans): 
	ncoils(n), radius(r), length(l), anglefunc(afunc), transfunc(tfunc)  {}
	
	void regularCoil(MixedSource& M, void* aparams = NULL, void* tparams = NULL) {
		buildEndpoints(aparams,tparams);
		buildSides(M);
		buildLineCaps(M);
	}
	
	const unsigned int ncoils;
	float radius;
	float length;
	mdouble (*anglefunc)(unsigned int, unsigned int, void*);
	vec3 (*transfunc)(int, int, int, int, int, void*);
	
	vec3 getEndp(unsigned int n, bool xside, bool yside, bool zside) { assert(n<ncoils); return endp[8*n + zside + 2*yside + 4*xside]; }
	
protected:
	
	std::vector<vec3> endp;
	
	void buildEndpoints(void* aparams, void* tparams) {
		mdouble x,y,phi;
		mdouble xmul,ymul,zmul;
		endp = std::vector<vec3>(8*ncoils);
		for(unsigned int i = 0; i<ncoils; i++) {
			for(int a1 = 0; a1<2; a1++) { // z
				for(int a2=0; a2<2; a2++) { // y
					for(int a3=0; a3<2; a3++) { // x
						xmul = (a3)?-1:1;
						ymul = (a2)?-1:1;
						zmul = (a1)?-1:1;
						phi = anglefunc(i+a3*ncoils+a2*2*ncoils+a1*4*ncoils,ncoils,aparams);
						x = radius*cos(phi);
						y = radius*sin(phi);
						endp[8*i + a1 + 2*a2 + 4*a3] = vec3(x*xmul,y*ymul,zmul*length*0.5) + transfunc(ncoils,i,a1,a2,a3,tparams);
					}
				}
			}
		}
	}
	
	void buildSides(MixedSource& M) {
		for(unsigned int i=0; i<ncoils; i++)
			for(unsigned int xside = 0; xside < 2; xside++)
				for(unsigned int yside = 0; yside < 2; yside++)
					M.addsource(new LineSource(getEndp(i,xside,yside,!yside),getEndp(i,xside,yside,yside),1.0/mdouble(ncoils)));
	}
	
	void buildArcCaps(MixedSource& M,unsigned int nseg=1) {
		
		for(unsigned int xside = 0; xside < 2; xside++) {
			for(unsigned int zside = 0; zside < 2; zside++) {
				
				for(unsigned int yside = 0; yside < 2; yside++) {
					for(unsigned int i=0; i<ncoils-1; i++) {
						if(yside == zside)
							M.arc(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)/mdouble(ncoils),nseg);
						else
							M.arc(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)/mdouble(ncoils),nseg);
					}
				}
				
				M.arc(getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),1.0,nseg);
			}
		}
		
	}
	
	void buildLineCaps(MixedSource& M) {
		
		for(unsigned int xside = 0; xside < 2; xside++) {
			for(unsigned int zside = 0; zside < 2; zside++) {
				
				for(unsigned int yside = 0; yside < 2; yside++) {
					for(unsigned int i=0; i<ncoils-1; i++) {
						if(yside == zside)
							M.addsource( new LineSource(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)/mdouble(ncoils)) );
						else
							M.addsource( new LineSource(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)/mdouble(ncoils)) );
					}
				}
				
				M.addsource( new LineSource( getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),1.0) );
			}
		}
		
	}
	
	void buildStraightCaps(MixedSource& M) {
		for(unsigned int xside = 0; xside < 2; xside++)
			for(unsigned int zside = 0; zside < 2; zside++)
				for(unsigned int i=0; i<ncoils; i++)
					M.addsource(new LineSource(getEndp(i,xside,zside,zside),getEndp(i,xside,!zside,zside),1.0/mdouble(ncoils)));
	}
	
	void buildMixedCaps(MixedSource& M, float rinner) {
		
		std::vector<vec3> innercircle[2][2][2];
		
		for(unsigned int xside = 0; xside < 2; xside++) {
			for(unsigned int zside = 0; zside < 2; zside++) {
				for(unsigned int i=0; i<ncoils; i++) {
					vec3 v1 = getEndp(i,xside,zside,zside);
					vec3 v2 = getEndp(i,xside,!zside,zside);
					if(fabs(v1[0]) >= rinner*radius) {
						M.addsource(new LineSource(v1,v2,1.0/mdouble(ncoils)));
					} else {
						float l = (1-sqrt(rinner*radius*rinner*radius-v1[0]*v1[0])/radius)*0.5;
						M.addsource(new LineSource(v1,v1*(1-l) + v2*l,1.0/mdouble(ncoils)));
						M.addsource(new LineSource(v1*l + v2*(1-l),v2,1.0/mdouble(ncoils)));
						innercircle[xside][0][zside].push_back(v1*(1-l) + v2*l);
						innercircle[xside][1][zside].push_back(v1*l + v2*(1-l));
					}
				}
			}
		}
		
		unsigned int imax = innercircle[0][0][0].size();
		vec3 dz[2];
		dz[0] = vec3(0,0,0.20*length);
		dz[1] = -dz[0];
		
		if(imax) {
			for(unsigned int xside = 0; xside < 2; xside++) {
				for(unsigned int zside = 0; zside < 2; zside++) {
					for(unsigned int yside = 0; yside < 2; yside++) {
						
						for(unsigned int i=0; i<imax-1; i++) {
							if(yside)
								M.arc(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i+1]+dz[zside],(i+1)/mdouble(ncoils),32);
							else
								M.arc(innercircle[xside][yside][zside][i+1]+dz[zside],innercircle[xside][yside][zside][i]+dz[zside],(i+1)/mdouble(ncoils),32);
						}
						
						for(unsigned int i=0; i<imax; i++) {
							if(!yside)
								M.addsource(new LineSource(innercircle[xside][yside][zside][i],innercircle[xside][yside][zside][i]+dz[zside],1/mdouble(ncoils)));
							else
								M.addsource(new LineSource(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i],1/mdouble(ncoils)));
						}
					}
					M.arc(innercircle[xside][true][zside][imax-1]+dz[zside],innercircle[xside][false][zside][imax-1]+dz[zside],mdouble(imax)/mdouble(ncoils),32);
				}
			}
		}
		
	}
	
	
};

#endif
