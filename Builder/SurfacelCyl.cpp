#include "SurfacelCyl.hh"
#include "ProgressBar.hh"
#include <cassert>

Surfacel_Protocol* Surfacel_Protocol::SP = new Surfacel_Protocol();

SurfacelSet::~SurfacelSet() {
	while(surfacels.size()) {
		surfacels.back()->release();
		surfacels.pop_back();
	}
}
	
vec3 SurfacelSet::fieldAt(const vec3& v) const {
	vec3 B;
	for(std::vector<ReactiveElement*>::const_iterator it = surfacels.begin(); it != surfacels.end(); it++)
		B += (*it)->fieldAt(v);
	return B;
}

void SurfacelSet::addSurfacel(ReactiveElement* e) {
	assert(nPhi==1);
	e->retain();
	surfacels.push_back(e);
	add_DF_group(e->nDF());
}

void SurfacelSet::visualize(bool top, mdouble scale) const {
	if(top) { vsr::startRecording(true); vsr::clearWindow(); }
	for(std::vector<ReactiveElement*>::const_iterator it = surfacels.begin(); it != surfacels.end(); it++)
		(*it)->visualize(false);
	if(top) vsr::stopRecording();
}

void SurfacelSet::calculateIncident(FieldSource* f) {
	if(verbose) printf("Calculating incident field on %i elements in %i groups...\n", n_subels(), n_subels()/nPhi);
	ProgressBar pb = ProgressBar(n_subels(), nPhi, verbose);
	incidentState = VarVec<mdouble>(nDF());
	for(unsigned int i = 0; i<n_subels(); i++) {
		pb.update(i);
		surfacels[i]->reactTo(f);
		for(unsigned int j = 0; j<surfacels[i]->nDF(); j++)
			incidentState[df_subindex(i,j)] = surfacels[i]->getState()[j];
	}
}

bool SurfacelSet::set_protocol(void* ip) {
	ReactiveUnitSet::set_protocol(ip);
	return ixn_ptcl==BField_Protocol::BFP || ixn_ptcl==Surfacel_Protocol::SP;
}

void SurfacelSet::queryInteraction() {
	if(ixn_ptcl == Surfacel_Protocol::SP)
		Surfacel_Protocol::SP->e = surfacels[ixn_el];
	else if(ixn_ptcl == BField_Protocol::BFP)
		BField_Protocol::BFP->B = surfacels[ixn_el]->fieldAt(BField_Protocol::BFP->x);
}

mvec SurfacelSet::subelReaction() {
	if(ixn_ptcl == Surfacel_Protocol::SP) {
		queryInteraction();
		return surfacels[ic_i]->interactionWith(Surfacel_Protocol::SP->e);
	} else if(ixn_ptcl == BField_Protocol::BFP) {
		assert(false);
	}
	return mvec(surfacels[ic_i]->nDF());
}

//----------------------------------------------

SurfacelCyl::~SurfacelCyl() {
	for(unsigned int i=0; i<segments.size(); i++)
		segments[i]->release();
}

void SurfacelCyl::append(ShieldSegment* s) {
	assert(s);
	s->retain();
	segments.push_back(s);
	for(unsigned int i=0; i<nPhi; i++) {
		surfacels.push_back(s->genElement(i));
		surfacels.back()->retain();
	}
	add_DF_group(surfacels.back()->nDF());
}

void SurfacelCyl::OptCone(unsigned int nZ0, unsigned int nZ1, vec2 s, vec2 e, PlanarElement* base, FieldEstimator2D* fes ) {
	
	//printf("Optimizing shield grid with %i fixed segments, %i varying segments\n",nZ0,nZ1);
	
	const unsigned int ngridpts = 2001;
	const unsigned int nZ = nZ0+nZ1;
	
	//cumulative field strength across length
	mdouble fstr[ngridpts];
	fstr[0]=0;
	vec2 dv = (e-s)*1./ngridpts;
	if(fes && nZ1) {
		for(unsigned int i=1; i<ngridpts; i++) {
			float l = (float(i)-0.5)/(ngridpts-1.0);
			//fstr[i] = fstr[i-1] + fes->estimateAt(s*(1-l)+e*l).mag(); 		// bunch up at high fields
			fstr[i] = fstr[i-1] + fabs(fes->derivAt(s*(1-l)+e*l,dv).dot(dv));	// bunch up at high field derivatives along sampling
		}
	} else {
		for(unsigned int i=0; i<ngridpts; i++)
			fstr[i] = i;
	}
	
	//printf("Average field %g\n",fstr[ngridpts-1]/(ngridpts-1));
	
	mdouble slope = 0;
	if(nZ1) {
		slope = fstr[ngridpts-1]*nZ0/((ngridpts-1.0)*nZ1); // add a constant slope for fixed partitions
		if(!(slope>1e-10)) slope=1e-10;	//< fix problems with extremely low fields
	}
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] += i*slope;
	
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] *= nZ/fstr[ngridpts-1];
	
	//interpolation to determine dividing lines, relative coordinates 0 to 1
	mdouble* ls = new mdouble[nZ+1];
	unsigned int n=1;
	int i=1;
	ls[0]=0;
	while(n<nZ) {
		if(fstr[i] >= n) {
			ls[n] = (mdouble(i)-(fstr[i]-n)/(fstr[i]-fstr[i-1]))/(ngridpts-1.0); // interpolated partition point
			n++;
			continue;
		}
		i++;	
	}
	ls[nZ] = 1.0;
	
	annulusSpec a;
	a.theta0 = 0;
	a.dTheta = 2*M_PI/mdouble(nPhi);
	
	base->retain();
	for(unsigned int z=0; z<nZ; z++) {
		a.start = s*(1-ls[z])+e*ls[z];
		a.end = s*(1-ls[z+1])+e*ls[z+1];
		append(new ShieldSegment(nPhi, base->reference(a)));
	}
	base->release();
	
	delete(ls);
}




