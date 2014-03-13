#include "SurfaceCurrentRS.hh"
#include "ProgressBar.hh"
#include <algorithm>

// SurfaceCurrent from interpolators of SurfaceCurrentRS
vec2 surfaceJ(vec2 v, void* p) {
	assert(p);
	return ((SurfaceCurrentRS*)p)->eval(v);
}

SurfaceCurrentRS::SurfaceCurrentRS(unsigned int nph, unsigned int nz): InterpolatingRS(nph), nZ(nz), Cref(NULL) {
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
	const SurfaceCurrentRS* S;	//< surface being integrated over
	vec2 ll0;					//< surface coordinate of active DF
	vec2 dlz;					//< offset of 1 unit in z coordinate
	vec2 dlp;					//< offset of 1 unit in phi coordinate
	int c_nz;					//< integrated element center nz
	int c_np;					//< integrated element center nphi
	bool pdir;					//< whether this is a response to current in the phi direction
	vec3 v;						//< position of interacting spot
	Matrix<2,3,mdouble> RM;		//< response transform to local field
};

mvec SRdA(vec2 l, void* params) {
	ResponseIntegParams& p = *(ResponseIntegParams*)params;
	
	vec2 lli = p.ll0 + vec2(l[0]/p.S->nZ, l[1]/p.S->nPhi);	//< surface coordinates of current source
	
	if(true) { //!p.c_nz || p.c_nz+1 == p.S->nZ) {
		// edge elements case --- don't try any shortcuts
		return mvec( p.RM * p.S->fieldAt_contrib_from(p.v, lli) );
	}
	
	// TODO make this agree with general method above
	
	vec3 xi = (*p.S->mySurface)(lli);		//< physical coordinates of current source
	vec3 r = xi - p.v;						//< r vector to integration pt.
	mdouble mr = r.mag();					//< distance between points
	if(mr < 1e-8) return mvec(2);
	
	vec3 dl = p.S->mySurface->deriv(lli,p.pdir);
	mdouble dw = p.S->mySurface->deriv(lli,!p.pdir).mag();
	mdouble interpl = p.S->Cref.basisFunc(l[0]-p.c_nz)*p.S->Cref.basisFunc(l[1]-p.c_np);
	
	return mvec( (p.RM * cross(dl,r)) * (interpl*dw/(4.*M_PI*mr*mr*mr)) );
}

vec2 SurfaceCurrentRS::subelReaction() {
	if(ixn_ptcl == BField_Protocol::BFP) {
		
		unsigned int el = ixn_df%(nZ*nPhi);
		vec2 sc = surf_coords(ic_i);
		
		// TODO proper self-interaction
		if(ic_i == el) return vec2(0,0);

		// set up integration info
		ResponseIntegParams p;
		p.S = this;
		p.v = (*mySurface)(sc);
		p.RM = sdefs[ic_i].rmat * mySurface->rotToLocal(sc);
		p.ll0 = surf_coords(0);
		p.c_nz = el/nPhi;
		p.c_np = el%nPhi;
		p.dlz = vec2(1./nZ,0);
		p.dlp = vec2(0,1./nPhi);
		p.pdir = ixn_df >= nZ*nPhi;
		
				
		int i_nz = ic_i/nPhi;
		int i_nphi = ic_i%nPhi;
		
		// save previous integration method
		Integration_Method im = myIntegrator.getMethod();
		
		//const unsigned int n_integ_domains = 3;						//< number of integration domains in each direction
		//const int integ_domains[n_integ_domains+1] = {-2,-1,1,2};	//< integration domain divisions
		// TODO allow 3x3 region with internal singularities
		
		const unsigned int n_integ_domains = 4;						//< number of integration domains in each direction
		const int integ_domains[n_integ_domains+1] = {-2,-1,0,1,2};	//< integration domain divisions


		// integrate over each of 16 interpolation regions
		vec2 vResp;
		for(unsigned int dmz = 0; dmz < n_integ_domains; dmz++) {
			int z0 = integ_domains[dmz];
			int z1 = integ_domains[dmz+1];
			if(p.c_nz + z1 < -1 || p.c_nz + z0 >= (int)nZ) continue; // skip past-edges domains
			for(unsigned int dmp = 0; dmp < n_integ_domains; dmp++) {
				int p0 = integ_domains[dmp];
				int p1 = integ_domains[dmp+1];
		
				int nz0 = p.c_nz + z0;
				int nz1 = p.c_nz + z1;
				int np0 = p.c_np + p0;
				int np1 = p.c_np + p1;
				
				// set integration method depending on whether there will be edge singularities
				if( (nz0 <= i_nz && i_nz <= nz1) && ( (i_nphi - np0 + nPhi)%nPhi <= (np1 - np0 + nPhi)%nPhi ) )
					myIntegrator.setMethod(INTEG_GSL_QAGS);
				else
					myIntegrator.setMethod(INTEG_GSL_QNG);
				
				mvec RB = myIntegrator.integrate2D(&SRdA, vec2(nz0<0 ? -0.5 : nz0, np0), vec2(nz1>=(int)nZ ? nZ-0.5 : nz1, np1), &p);
				
				//std::cout << vec2(nz0,np0) << vec2(nz1,np1) << " " << RB << std::endl;
				
				vResp += vec2(RB[0],RB[1]);
			}
		}
		
		// restore previous integration method
		myIntegrator.setMethod(im);
		
		//std::cout << vec2(i_nz,i_nPhi) << " to " << vec2(p.c_nz,p.c_nphi) << p.pdir << ":\t" << vResp / (nZ * nPhi) << std::endl;
		
		return vResp / (nZ * nPhi);	//< note adjust for integrating over [0,nZ]*[0,nPhi] insead of [0,1]^2
		
	}
	assert(false);
	return vec2();
}

mdouble SurfaceCurrentRS::nextInteractionTerm(unsigned int& i, unsigned int& j) {
	// get previous cached term
	mdouble v = ic_v[ic_di];
	i = ic_i + nZ*nPhi*ic_di;
	j = ixn_df;
	
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
