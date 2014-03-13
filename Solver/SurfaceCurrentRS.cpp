#include "SurfaceCurrentRS.hh"
#include "ProgressBar.hh"
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
	
	// split up surface phi integrals by default
	dflt_integrator_ndivs_y = nPhi>8? nPhi/8 : 1;
	// use adaptive integration by default
	myIntegrator.setMethod(INTEG_GSL_QAG);
}

mvec SurfaceCurrentRS::getReactionTo(ReactiveSet* R, unsigned int phi) {
	mvec v(nDF()/nPhi);
	for(ixn_el = phi; ixn_el < nDF()/2; ixn_el+=nPhi) {
		vec2 vs = subelReaction(R);
		v[ixn_el/nPhi] = vs[0];
		v[ixn_el/nPhi + nZ] = vs[1];
	}
	return v;
}

vec2 SurfaceCurrentRS::subelReaction(ReactiveSet* R) {
	vec2 sc = surf_coords(ixn_el);
	BField_Protocol::BFP->x = (*mySurface)(sc);
	Matrix<2,3,mdouble> RM = sdefs[ixn_el].rmat * mySurface->rotToLocal(sc);
	BField_Protocol::BFP->M = &RM;
	BField_Protocol::BFP->caller = this;
	if(!R->queryInteraction(BField_Protocol::BFP)) { assert(false); return vec2(0,0); }
	return BField_Protocol::BFP->MB;
}
	
bool SurfaceCurrentRS::queryInteraction(void* ip) {

	if(ip != BField_Protocol::BFP) return false;
	
	BField_Protocol::BFP->B = vec3(0,0,0);
	BField_Protocol::BFP->MB = vec2(0,0);
	
	// z and phi numbers of active and responding element
	unsigned int el = ixn_df % (nZ*nPhi);
	if(BField_Protocol::BFP->caller == this && el == ixn_el) return true; // TODO proper self-interaction
	int c_nz = el/nPhi;
	int c_np = el%nPhi;
	int i_nz = ixn_el/nPhi;
	int i_nphi = ixn_el%nPhi;
	
	// How to slice up integration range; TODO allow 3x3 region with internal singularities
	//const unsigned int n_integ_domains = 3;					//< number of integration domains in each direction
	//const int integ_domains[n_integ_domains+1] = {-2,-1,1,2};	//< integration domain divisions
	const unsigned int n_integ_domains = 4;						//< number of integration domains in each direction
	const int integ_domains[n_integ_domains+1] = {-2,-1,0,1,2};	//< integration domain divisions

	// save previous integration method
	Integration_Method im = myIntegrator.getMethod();
	
	// integrate over each region
	vec2 vResp;
	for(unsigned int dmz = 0; dmz < n_integ_domains; dmz++) {
		int z0 = integ_domains[dmz];
		int z1 = integ_domains[dmz+1];
		if(c_nz + z1 < 0 || c_nz + z0 >= (int)nZ) continue; // skip past-edges domains
		for(unsigned int dmp = 0; dmp < n_integ_domains; dmp++) {
			int p0 = integ_domains[dmp];
			int p1 = integ_domains[dmp+1];
	
			int nz0 = c_nz + z0;
			int nz1 = c_nz + z1;
			int np0 = c_np + p0;
			int np1 = c_np + p1;
			
			// set integration method depending on whether there will be edge singularities
			if( BField_Protocol::BFP->caller == this && (nz0 <= i_nz && i_nz <= nz1) && ( (i_nphi - np0 + nPhi)%nPhi <= (np1 - np0 + nPhi)%nPhi ) )
				myIntegrator.setMethod(INTEG_GSL_QAGS);
			else
				myIntegrator.setMethod(INTEG_GSL_QAG);
			
			vec2 ll(double(nz0)/double(nZ), double(np0)/double(nPhi));
			vec2 ur(double(nz1)/double(nZ), double(np1)/double(nPhi));
						
			if(BField_Protocol::BFP->M) {
				BField_Protocol::BFP->MB += fieldAtWithTransform(BField_Protocol::BFP->x, *BField_Protocol::BFP->M, ll, ur, 1, 1);
			} else {
				BField_Protocol::BFP->B += fieldAt(BField_Protocol::BFP->x, ll, ur, 1, 1);
			}
		}
	}
	
	// restore previous integration method
	myIntegrator.setMethod(im);
	
	return true;
}


//
// Response to incident magnetic fields
//

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

struct AverageFieldIntegParams {
	const FieldSource* f;			//< field source
	const SurfaceCurrentRS* S;		//< this surface being integrated over
	Matrix<2,3,mdouble> rmat;		//< surface response to local field
};

mvec FieldResponsedA(vec2 l, void* params) {
	AverageFieldIntegParams* p = (AverageFieldIntegParams*)params;
	vec2 r = p->rmat * p->S->mySurface->rotToLocal(l) * p->f->fieldAt((*p->S->mySurface)(l));
	return mvec(r) * p->S->mySurface->dA(l);
}

void SurfaceCurrentRS::calculateIncident(const FieldSource& f) {
	assert(sdefs.size() == nZ*nPhi);
	assert(mySurface);
	
	incidentState = mvec(nDF());
	
	// save previous integration method
	Integration_Method im = myIntegrator.getMethod();
	myIntegrator.setMethod(INTEG_GSL_QAG);
	
	printf("Calculating incident field over %i elements...\n",nZ*nPhi);
	ProgressBar pb = ProgressBar(nZ, 1, true);
		
	for(unsigned int zn = 0; zn < nZ; zn++) {
		pb.update(zn);
		for(unsigned int pn = 0; pn < nPhi; pn++) {
			
			unsigned int el = zn*nPhi + pn;
	
			AverageFieldIntegParams AFIP;
			AFIP.f = &f;
			AFIP.S = this;
			AFIP.rmat = sdefs[el].rmat;
		
			vec2 ll(double(zn)/double(nZ), double(pn)/double(nPhi));
			vec2 ur(double(zn+1)/double(nZ), double(pn+1)/double(nPhi));
			
			mvec r = myIntegrator.integrate2D(&FieldResponsedA, ll, ur, &AFIP);
			double A = mySurface->area(ll,ur);
			incidentState[el] = r[0]/A;
			incidentState[el+nZ*nPhi] = r[1]/A;
			
			// value at center:
			//vec2 l = (ll+ur)*0.5;
			//vec2 r = AFIP.rmat * AFIP.S->mySurface->rotToLocal(l) * AFIP.f->fieldAt((*AFIP.S->mySurface)(l));
			//incidentState[el] = r[0];
			//incidentState[el+nZ*nPhi] = r[1];
		}
		
	}
	
	// restore previous integration method
	myIntegrator.setMethod(im);
			
	_setDF(incidentState);
}
