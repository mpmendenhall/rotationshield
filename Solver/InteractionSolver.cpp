#include "InteractionSolver.hh"
#include "ProgressBar.hh"

InteractionSolver::InteractionSolver(unsigned int gs): ndf(0), groupSize(gs), groupsComplete(true), fInc(NULL),
surfacels(std::vector<ReactiveElement*>()), groupIndex(std::vector<unsigned int>()), verbose(false) {
	groupIndex.push_back(0);
}

InteractionSolver::~InteractionSolver() { 
	while(N()) {
		surfacels.back()->release(); 
		surfacels.pop_back();
	}
	if(fInc) 
		fInc->release();
}

void InteractionSolver::addSurfacel(ReactiveElement* e) {
	e->retain();
	surfacels.push_back(e);
	ndf += e->nDF();
	groupsComplete = !(N()%groupSize);
	if(groupsComplete)
		groupIndex.push_back(ndf);
}

void InteractionSolver::calculateIncident(FieldSource* f) {
	assert(groupsComplete);
	if(verbose) printf("Calculating incident field on %i elements in %i groups...\n",N(),N()/groupSize);
	ProgressBar pb = ProgressBar(N(),groupSize,verbose);
	if(fInc) fInc->release();
	fInc = f;
	fInc->retain();
	incidentState = VarVec<mdouble>(ndf);
	for(unsigned int i = 0; i<N(); i++) {
		pb.update(i);
		surfacels[i]->reactTo(fInc);
		for(unsigned int j = 0; j<surfacels[i]->nDF(); j++)
			incidentState[index(i,j)] = surfacels[i]->getState()[j];
	}
}	

vec3 InteractionSolver::fieldAt(const vec3& v) const { 
	vec3 r = vec3();
	for(unsigned int i=0; i<N(); i++)
		r += surfacels[i]->fieldAt(v);
	return r;
}

vec3 InteractionSolver::partialField(const vec3& v, unsigned int n1, int nf) const {
	unsigned int n2 = (nf<0)?N():abs(nf);
	assert(n1<=n2 && n2<=N());
	vec3 r = vec3();
	for(unsigned int i=n1; i<n2; i++)
		r += surfacels[i]->fieldAt(v);
	return r;
}

void InteractionSolver::visualize(bool top, mdouble scaling) const {
	if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); }
	for(unsigned int i=0; i<N(); i++)
		surfacels[i]->visualize(false);
	if(top) vsr::Visr::W->stopRecording();
}

mvec InteractionSolver::getFinalState(unsigned int i) const {
	assert(i<N());
	assert(finalState.size() == ndf);
	mvec s = mvec(surfacels[i]->nDF());
	for(unsigned int j=0; j<surfacels[i]->nDF(); j++)
		s[j] = finalState[index(i,j)];
	return s;
}

mmat InteractionSolver::getInteraction(int x0, int x1) const {
	if(x0 == x1) {
		mmat mInt = mmat(surfacels[x0]->nDF(),surfacels[x0]->nDF());
		for(unsigned int i=0; i<surfacels[x0]->nDF(); i++)
			mInt(i,i) = 1.0;
		return mInt;
	}
	return -surfacels[x0]->interactionWith(surfacels[x1]);
}
