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
	std::vector<InterpolationHelper*> iv = InterplDF.getSubHelpers(1);
	for(std::vector<InterpolationHelper*>::iterator it = iv.begin(); it != iv.end(); it++)
		(*it)->getInterpolator()->setSymmetricOffset();
	iv = InterplDF.getSubHelpers(2);
	for(std::vector<InterpolationHelper*>::iterator it = iv.begin(); it != iv.end(); it++)
		(*it)->getInterpolator()->setSymmetricOffset();
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
	assert(mySurface);
	for(unsigned int el = 0; el < nZ*nPhi; el++)
		vis_coords(surf_coords(el));
}

vec2 SurfaceCurrentRS::eval(const vec2& p) const {
	const unsigned int zero = 0;
	const unsigned int one = 1;
	mdouble jx = InterplDF.getSubHelper(1,&zero).eval(&p[0]);
	mdouble jy = InterplDF.getSubHelper(1,&one).eval(&p[0]);
	return vec2(jx,jy);
}

void SurfaceCurrentRS::calculateIncident(const FieldSource& f) {
	assert(sdefs.size() == nZ*nPhi);
	assert(mySurface);
	
	incidentState = mvec(nDF());
	
	for(unsigned int el = 0; el < nZ*nPhi; el++) {
		vec2 l = surf_coords(el);
		vec3 x = (*mySurface)(l);
		vec3 Bl = mySurface->rotToLocal(l)*f.fieldAt(x);
		vec2 r = sdefs[el].responseToField(Bl);
	
		incidentState[el] = r[0];
		incidentState[el+nZ*nPhi] = r[1];
	}
	
	_setDF(incidentState);
}
