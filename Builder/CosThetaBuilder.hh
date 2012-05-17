#ifndef COSTHETABUILDER_HH
#define COSTHETABUILDER_HH 1


#include "MixedSource.hh"
#include "MiscUtils.hh"
#include <vector>

/// Cos theta coil with fourier-series of wire shifts from normal
mdouble shiftPositioner(unsigned int i, unsigned int ncoils, void* params);
/// Cos theta coil positioner using "k" parameter
mdouble alarconKPositioner(unsigned int i, unsigned int ncoils, void* params);

vec3 simpletrans(int N, int i, int a1, int a2, int a3, void* params);

struct doubleTrans {
	vec3 trans1;
	vec3 trans2;
};

/// different translation for +z,-z ends
vec3 fancytrans(int N, int i, int a1, int a2, int a3, void* params);

class CosThetaBuilder {
public:
	CosThetaBuilder(unsigned int n, mdouble r, mdouble l,
					mdouble (*afunc)(unsigned int, unsigned int, void*) = &shiftPositioner,
					vec3 (*tfunc)(int, int, int, int, int, void*) = &simpletrans): 
	ncoils(n), radius(r), length(l), anglefunc(afunc), transfunc(tfunc)  {}
	
	void regularCoil(MixedSource& M, void* aparams = NULL, void* tparams = NULL);
	
	const unsigned int ncoils;
	float radius;
	float length;
	
	mdouble (*anglefunc)(unsigned int, unsigned int, void*);
	vec3 (*transfunc)(int, int, int, int, int, void*);
	vec3 getEndp(unsigned int n, bool xside, bool yside, bool zside);
	
protected:
	
	std::vector<vec3> endp;
	
	void buildEndpoints(void* aparams, void* tparams);	
	void buildSides(MixedSource& M);	
	void buildArcCaps(MixedSource& M,unsigned int nseg=1);	
	void buildLineCaps(MixedSource& M);	
	void buildStraightCaps(MixedSource& M);	
	void buildMixedCaps(MixedSource& M, float rinner);	
};

#endif
