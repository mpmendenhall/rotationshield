#include "SurfaceCurrentRS.hh"
#include <algorithm>

// SurfaceCurrent from interpolators of SurfaceCurrentRS
vec2 surfaceJ(vec2 v, void* p) {
	assert(p);
	SurfaceCurrentRS& S = *(SurfaceCurrentRS*)p;
	
	return vec2();
}

SurfaceCurrentRS::SurfaceCurrentRS(unsigned int nph, unsigned int nz): InterpolatingRS(nph), nZ(nz) {
	// surface current function
	sj = &surfaceJ;
	sjparams = this;
	
	// set up DF interpolator
	unsigned int ndm[] = {2,nZ,nPhi};
	InterplDF.setupDataGrid(ndm, ndm+3);
}

bool SurfaceCurrentRS::set_protocol(void* ip) {
	InterpolatingRS::set_protocol(ip);
	return ixn_ptcl==BField_Protocol::BFP;
}

void SurfaceCurrentRS::queryInteraction() {
	if(ixn_ptcl == BField_Protocol::BFP)
		BField_Protocol::BFP->B = fieldAt(BField_Protocol::BFP->x);
}

void SurfaceCurrentRS::setSurfaceResponse(SurfaceI_Response r) {
	sdefs.resize(nZ*nPhi);
	std::fill(sdefs.begin(), sdefs.end(), r);
}

