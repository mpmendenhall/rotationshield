#include "SurfaceCurrentRS.hh"
#include <algorithm>

// SurfaceCurrent from interpolators of SurfaceCurrentRS
vec2 surfaceJ(vec2 v, void* p) {
	assert(p);
	return ((SurfaceCurrentRS*)p)->eval(v);
}

SurfaceCurrentRS::SurfaceCurrentRS(unsigned int nph, unsigned int nz): InterpolatingRS(nph), nZ(nz) {
	// surface current function
	sj = &surfaceJ;
	sjparams = this;
	
	// set up DF interpolator
	unsigned int ndm[] = {2,nZ,nPhi};
	InterplDF.setupDataGrid(ndm, ndm+3);
	InterplDF.setBoundaryCondition(BC_CYCLIC,2);
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

void SurfaceCurrentRS::_visualize() const {
	SurfaceCurrentSource::_visualize();
	// TODO show element vectors
}

vec2 SurfaceCurrentRS::eval(const vec2& p) const {
	mdouble jx = InterplDF.getSubHelpers()[0]->eval(&p[0]);
	mdouble jy = InterplDF.getSubHelpers()[1]->eval(&p[0]);
	return vec2(jx,jy);
}

