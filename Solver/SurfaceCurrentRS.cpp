#include "SurfaceCurrentRS.hh"
#include "ProgressBar.hh"
#include "VisSurface.hh"
#include <algorithm>

void SurfaceCurrentRS::element_surface_range(unsigned int i, vec2& ll, vec2& ur) const {
	unsigned int zn = i/nPhi;
	unsigned int pn = i%nPhi;
	ll = vec2( double(zn)/nZ, double(pn)/nPhi );
	ur = vec2( double(zn+1)/nZ, double(pn+1)/nPhi );
}


// SurfaceCurrent from interpolators of SurfaceCurrentRS
vec2 surfaceJ(vec2 v, void* p) {
	assert(p);
	return ((SurfaceCurrentRS*)p)->eval(v);
}

SurfaceCurrentRS::SurfaceCurrentRS(unsigned int nph, unsigned int nz): InterpolatingRS2D(nph) {
	// surface current function
	sj = &surfaceJ;
	sjparams = this;
	
	// set up DF interpolator
	make_grids(nz,2);
				
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
	unsigned int el = ixn_df % (nZ*nPhi);	//< active element
	if(BField_Protocol::BFP->caller == this && el == ixn_el) return true; // TODO proper self-interaction
	int c_nz = el/nPhi;			//< center z of active element
	int c_np = el%nPhi;			//< center phi of active element
	
	int i_nz = ixn_el/nPhi;		//< responding element z
	int i_nphi = ixn_el%nPhi;	//< responding element phi
	
	
	// How to slice up integration range; TODO allow 3x3 region with internal singularities
	const unsigned int n_integ_domains = 3;					//< number of integration domains in each direction
	const int integ_domains[n_integ_domains+1] = {-2,-1,1,2};	//< integration domain divisions
	//const unsigned int n_integ_domains = 4;						//< number of integration domains in each direction
	//const int integ_domains[n_integ_domains+1] = {-2,-1,0,1,2};	//< integration domain divisions

	// save previous integration method
	Integration_Method im = myIntegrator.getMethod();
	
	// integrate over each region
	for(unsigned int dmz = 0; dmz < n_integ_domains; dmz++) {
		for(unsigned int dmp = 0; dmp < n_integ_domains; dmp++) {

			int nz0 = c_nz + integ_domains[dmz];
			int nz1 = c_nz + integ_domains[dmz+1];
			int np0 = c_np + integ_domains[dmp];
			int np1 = c_np + integ_domains[dmp+1];
			if(nz1 <= 0 || nz0 >= (int)nZ) continue; // skip past-edges domains; note, we are allowed to wrap around in phi
			
			vec2 ll(double(nz0)/double(nZ), double(np0)/double(nPhi));
			vec2 ur(double(nz1)/double(nZ), double(np1)/double(nPhi));
			if(ll[0]<0) ll[0]=0;
			if(ur[0]>1) ur[0]=1;
			
			// set integration method depending on whether there will be edge singularities (target point inside integration range)
			if( BField_Protocol::BFP->caller == this
					&& (nz0 <= i_nz && i_nz <= nz1)
					&& ( (i_nphi - np0 + nPhi)%nPhi <= (np1 - np0 + nPhi)%nPhi ) )
				myIntegrator.setMethod(INTEG_GSL_QAGS);
			else
				myIntegrator.setMethod(INTEG_GSL_QAG);
												
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
	vsr::setColor(0,0,0);
	for(unsigned int el = 0; el < nZ*nPhi; el++)
		vis_coords(surf_coords(el));
	vis_i_vectors();
}

void SurfaceCurrentRS::vis_i_vectors(double s, double mx) const {
	for(unsigned int el = 0; el < nZ*nPhi; el++) {
		vec2 l = surf_coords(el);
		vec3 o = (*mySurface)(l);
		
		vec2 ll,ur;
		element_surface_range(el, ll, ur);
		vec3 j = netCurrent(ll, ur, 1, 1)*s/mySurface->area(ll,ur);
		double m = j.mag();
		j *= mx*atan(s*m/mx)/(m*0.5*M_PI);
		
		vsr::setColor(0,0,1);
		vsr::line(o-j*0.5, o);
		vsr::setColor(1,0,0);
		vsr::line(o, o+j*0.5);
		
		/*
		vec3 dx = mySurface->deriv(l,0).normalized();
		vec3 dy = mySurface->deriv(l,1).normalized();
		vec2 j = eval(l)*s;
		vsr::line(o, o+dx*j[0]+dy*j[1]);
		*/
	}
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
			
	_setDFv(incidentState);
}
