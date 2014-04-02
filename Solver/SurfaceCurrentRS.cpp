/* 
 * SurfaceCurrentRS.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

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
	return ((SurfaceCurrentRS*)p)->eval_J(v);
}

SurfaceCurrentRS::SurfaceCurrentRS(SurfaceGeometry* SG, unsigned int nph, unsigned int nz, const std::string& nm):
SurfaceCurrentSource(SG,nm), InterpolatingRS2D(nph), point_ixn(true) {
	
	assert(mySurface);
	assert(nPhi == 1 || mySurface->isClosed(1));
	
	// surface current function
	sj = &surfaceJ;
	sjparams = this;
	
	// set up DF interpolator
	make_grids(nz,2);
	if(mySurface->isClosed(0))
		for(unsigned int i=0; i<nDFi; i++)
			G[i]->bc[0] = IB_CYCLIC;
	
	// split up surface phi integrals by default
	mySurface->dflt_integrator_ndivs_y = nPhi>4? nPhi/4 : 1;
}


void SurfaceCurrentRS::set_current_loop(unsigned int z, double i, bool phidir) {
	if(phidir) {
		assert(z<nZ);
		for(unsigned int p=0; p<nPhi; p++)
			setDF(nZ*nPhi + z*nPhi + p, i);
	} else {
		assert(z<nPhi);
		for(unsigned int p=0; p<nZ; p++)
			setDF(p*nPhi + z, i);
	}
}

mvec SurfaceCurrentRS::getReactionTo(ReactiveSet* R, unsigned int phi) {
	mvec v(nDF()/nPhi);
	for(ixn_el = phi; ixn_el < nDF()/nDFi; ixn_el+=nPhi) {
		mvec vs = subelReaction(R);
		if(vs.size() != nDFi) { assert(false); continue; }
		for(unsigned int dfi = 0; dfi < nDFi; dfi++)
			v[ixn_el/nPhi + nZ*dfi] = vs[dfi];
	}
	return v;
}

/// parameters for averaging field over a region
struct AverageFieldIntegParams {
	ReactiveSet* RS;			//< field source
	const SurfaceCurrentRS* S;	//< this surface being integrated over
	Matrix<2,3,double> rmat2;	//< 2-component response matrix to local field
};

/// field response integration function via protocol
mvec fieldResponse_IntegF(mvec v, void* params) {
	AverageFieldIntegParams* p = (AverageFieldIntegParams*)params;
	vec2 l(v[0],v[1]);
	BField_Protocol::BFP->x = (*p->S->mySurface)(l);
	Matrix<2,3,double> RM2 = p->rmat2 * p->S->mySurface->rotToLocal(l);
	BField_Protocol::BFP->M2 = &RM2;
	BField_Protocol::BFP->M2B = vec2(0,0);
	if(!p->RS->queryInteraction(BField_Protocol::BFP)) { assert(false); return mvec(2); }
	return mvec(BField_Protocol::BFP->M2B * p->S->mySurface->dA(l));
}

mvec SurfaceCurrentRS::subelReaction(ReactiveSet* R) {
	if(R->interactionMode() || point_ixn) {
		vec2 sc = surf_coords(ixn_el);
		BField_Protocol::BFP->x = (*mySurface)(sc);
		Matrix<2,3,double> RM2 = sdefs[ixn_el].rmat2 * mySurface->rotToLocal(sc);
		BField_Protocol::BFP->M2 = &RM2;
		BField_Protocol::BFP->M2B = vec2(0,0);
		BField_Protocol::BFP->M3 = NULL;
		BField_Protocol::BFP->caller = RS_UID;
		if(!R->queryInteraction(BField_Protocol::BFP)) { assert(false); return mvec(); }
		return mvec(BField_Protocol::BFP->M2B);
	} else {
		BField_Protocol::BFP->M3 = NULL;
		BField_Protocol::BFP->caller = RS_UID;
		
		AverageFieldIntegParams AFIP;
		AFIP.RS = R;
		AFIP.S = this;
		AFIP.rmat2 = sdefs[ixn_el].rmat2;
		
		vec2 ll,ur;
		element_surface_range(ixn_el, ll, ur);
		mvec r = mySurface->myIntegratorND.integrate(&fieldResponse_IntegF, 2, ll, ur, &AFIP);
		double A = mySurface->area(ll,ur);
		return r/A;
	}
}
	
bool SurfaceCurrentRS::queryInteraction(void* ip) {

	if(ip != BField_Protocol::BFP) return false;
	
	// case for all DF set; respond with magnetic field
	if(ixn_df >= nDF()) {
		if(BField_Protocol::BFP->M2) {
			BField_Protocol::BFP->M2B += fieldAtWithTransform2(BField_Protocol::BFP->x, *BField_Protocol::BFP->M2);
		} else if(BField_Protocol::BFP->M3) {
			BField_Protocol::BFP->B += fieldAtWithTransform3(BField_Protocol::BFP->x, *BField_Protocol::BFP->M3);
		} else {
			BField_Protocol::BFP->B += fieldAt(BField_Protocol::BFP->x);
		}
		return true;
	}
	
	// z and phi numbers of active and responding element
	unsigned int el = ixn_df % (nZ*nPhi);	//< active element
	int c_nz = el/nPhi;			//< center z of active element
	int c_np = el%nPhi;			//< center phi of active element
	
	int i_nz = ixn_el/nPhi;		//< responding element z
	int i_nphi = ixn_el%nPhi;	//< responding element phi
	vec2 ixn_center = surf_coords(ixn_el);	//< responding element coordinate center
	
	// whether this is the interaction of an element with itself
	//bool self_ixn = BField_Protocol::BFP->caller == RS_UID && el == ixn_el;

	// distance between active and responding element
	unsigned int delta_phi = (i_nphi-c_np + 2*nPhi)%nPhi;
	if(delta_phi >= nPhi/2) delta_phi -= nPhi;
	int delta_z = i_nz-c_nz;
	if(mySurface->isClosed(0)) {
		delta_z = (delta_z + 2*nZ)%nZ;
		if(delta_z >= int(nZ)/2) delta_z -= nZ;
	}
	
	bool self_intersection = (BField_Protocol::BFP->caller == RS_UID && delta_z == 0 && delta_phi == 0);
	
	// How to slice up integration range
	std::vector<int> integ_domains;
	const int idomains_9[] = {-2,-1,1,2};
	integ_domains.insert(integ_domains.end(), idomains_9, idomains_9+4);
	
	mySurface->polar_r0 = 0;
	unsigned int nerrx = mySurface->myIntegrator2D.reset_errcount();
	unsigned int nerry = mySurface->myIntegrator2D.reset_y_errcount();
	
	// integrate over each region
	for(unsigned int dmz = 0; dmz < integ_domains.size()-1; dmz++) {
		for(unsigned int dmp = 0; dmp < integ_domains.size()-1; dmp++) {

			int nz0 = c_nz + integ_domains[dmz];
			int nz1 = c_nz + integ_domains[dmz+1];
			int np0 = c_np + integ_domains[dmp];
			int np1 = c_np + integ_domains[dmp+1];
			// skip past-edges domains; note, we are allowed to wrap around in phi
			if(!mySurface->isClosed(0) && (nz1 < 0 || nz0 >= (int)nZ)) continue;
			
			bool atcenter = integ_domains[dmz] <= 0 && 0 <= integ_domains[dmz+1] && integ_domains[dmp] <= 0 && 0 <= integ_domains[dmp+1];
			if(self_intersection && atcenter) {
				mySurface->myIntegrator2D.setMethod(INTEG_GSL_CQUAD);
				mySurface->polar_integral_center = &ixn_center;
				if(np0 > i_nphi) { np0 -= nPhi; np1 -= nPhi; }						// adjust phi definition to align with polar center
				if(mySurface->isClosed(0) && nz0 > i_nz) { nz0 -= nZ; nz1 -= nZ; }	// adjust z definition to align with polar center
			}
			
			vec2 ll((nz0+0.5)/double(nZ), (np0+0.5)/double(nPhi));
			vec2 ur((nz1+0.5)/double(nZ), (np1+0.5)/double(nPhi));
			if(!mySurface->isClosed(0)) {
				if(ll[0]<0) ll[0]=0;
				if(ur[0]>1) ur[0]=1;
			}
			
			if(!(self_intersection && atcenter)) {
				// select integration method depending on near/far difference
				double mn,mx;
				mySurface->proximity(BField_Protocol::BFP->x, ll, ur, mn, mx);
				if(mx > 4*mn) mySurface->myIntegrator2D.setMethod(INTEG_GSL_CQUAD);
				else mySurface->myIntegrator2D.setMethod(INTEG_GSL_QAG);
			}
						
			if(BField_Protocol::BFP->M2) {
				BField_Protocol::BFP->M2B += fieldAtWithTransform2(BField_Protocol::BFP->x, *BField_Protocol::BFP->M2, ll, ur, 1, 1);
			} else if(BField_Protocol::BFP->M3) {
				BField_Protocol::BFP->B += fieldAtWithTransform3(BField_Protocol::BFP->x, *BField_Protocol::BFP->M3, ll, ur, 1, 1);
			} else {
				BField_Protocol::BFP->B += fieldAt(BField_Protocol::BFP->x, ll, ur, 1, 1);
			}
				
			nerrx = mySurface->myIntegrator2D.reset_errcount();
			nerry = mySurface->myIntegrator2D.reset_y_errcount();
			if(nerrx || nerry) {
				std::cout << "range:" << ll << ur << " active el:" << el << " " << ixn_center
					<< " P" << bool(mySurface->polar_integral_center) << " method " << mySurface->myIntegrator2D.getMethod()
					<< " from " << (*mySurface)(surf_coords(ixn_df)) << " to " << BField_Protocol::BFP->x
					<< " asking:" << BField_Protocol::BFP->caller << " responding:" << RS_UID
					<< " errs (" << nerrx << "," << nerry
					<< ") z: " << nz0 << "/" << i_nz << "/" << nz1 << " phi: " << np0 << "/" << i_nphi << "/" << np1 << std::endl;
			}
			
			mySurface->polar_integral_center = NULL;
		}
	}
		
	mySurface->myIntegrator2D.setMethod(INTEG_GSL_QAG); //< reset default method
	
	return true;
}

void SurfaceCurrentRS::setSurfaceResponse(SurfaceI_Response r) {
	sdefs.resize(nZ*nPhi);
	std::fill(sdefs.begin(), sdefs.end(), r);
}

void SurfaceCurrentRS::_visualize() const {
	SurfaceCurrentSource::_visualize();
	assert(mySurface);
	vsr::setColor(0,0,0);
	for(unsigned int el = 0; el < nZ*nPhi; el++) {
		vec2 l = surf_coords(el);
		vis_coords(l);
	}
	vis_i_vectors();
}

void SurfaceCurrentRS::vis_i_vectors(double s, double mx) const {
	if(!mySurface || !sj) return;
	
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
	}
}

