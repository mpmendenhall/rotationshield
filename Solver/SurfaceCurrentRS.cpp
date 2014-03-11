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
		
	set_protocol(BField_Protocol::BFP);
}

bool SurfaceCurrentRS::set_protocol(void* ip) {
	InterpolatingRS::set_protocol(ip);
	return ixn_ptcl==BField_Protocol::BFP;
}

void SurfaceCurrentRS::queryInteraction() {
	if(ixn_ptcl == BField_Protocol::BFP) {
		unsigned int el = ixn_df%(nZ*nPhi);
		vec2 l = surf_coords(el);
		vec2 dl(2.1/nZ, 2.1/nPhi);
		BField_Protocol::BFP->B = fieldAt(BField_Protocol::BFP->x, l-dl, l+dl);
	}
}


struct ResponseIntegParams {
	const SurfaceSource* S;
	vec3 v;					//< position in local field
	Matrix<2,3,mdouble> RM;	//< response transform to local field
	
};

mvec SRdA(mdouble x, mdouble y, void* params) {
	ResponseIntegParams& p = *(ResponseIntegParams*)params;
	return mvec(p.RM * p.S->fieldAt_contrib_from(p.v,x,y));
}

mdouble clamp1(mdouble m) { return m<0? 0 : m<1? m : 1; }

vec2 SurfaceCurrentRS::subelReaction() {
	if(ixn_ptcl == BField_Protocol::BFP) {
		vec2 l = surf_coords(ic_i);	//< position on this surface
	
		// warn integrator about singular point
		myIntegrator.xysingularities.clear();
		myIntegrator.xysingularities.push_back(l);
		
		ResponseIntegParams p;
		p.S = this;
		p.v = (*mySurface)(l);
		p.RM = sdefs[ic_i].rmat * mySurface->rotToLocal(l);
		
		vec2 dl(2.01/nZ, 2.01/nPhi);
		vec2 ll = l-dl;
		vec2 ur = l+dl;
		mvec RB = myIntegrator.integrate(&SRdA, clamp1(ll[0]), clamp1(ur[0]), clamp1(ll[1]), clamp1(ur[1]), &p);
		
		return vec2(RB[0],RB[1]);
	}
	assert(false);
	return vec2();
}

mdouble SurfaceCurrentRS::nextInteractionTerm(unsigned int& i, unsigned int& j) {
	// get previous cached term
	mdouble v = ic_v[ic_di];
	i = ic_i + nZ*nPhi*ic_di;
	j = ixn_df;
	
	// self-interaction
	if(i==j) v = 0;
	
	// step to next term
	ic_di = (ic_di+1)%2;
	if(!ic_di) {
		setInteractionDF((ixn_df+1)%nDF());
		if(!ixn_df) {
			ic_i = (ic_i+nPhi)%(nZ*nPhi);
			if(ic_i < nPhi) ic_i = (ic_i+1)%nPhi;
		}
		ic_v = subelReaction();
	}
	
	return v;
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
