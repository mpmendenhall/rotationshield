#include "ReactiveSet.hh"
#include <iostream>

void ReactiveSet::setInteractionDF(unsigned int DF, double v) {

	//std::cout << this << " Setting Interaction DF " << DF << " of " << nDF() << " to " << v << std::endl;

	assert(DF < finalState.size() && ixn_df < finalState.size());
	finalState[ixn_df] = 0;
	finalState[DF] = v;
	_setDF(ixn_df,0);
	_setDF(DF,v);
	ixn_df = DF;
}

void ReactiveSet::_setDFv(const mvec& v) {
	assert(v.size() == nDF());
	for(unsigned int DF = 0; DF < nDF(); DF++)
		_setDF(DF,v[DF]);
}


//-----------------------------------------

mvec ReactiveUnitSet::getReactionTo(ReactiveSet* R, unsigned int phi) {
	mvec v(nDF()/nPhi);
	for(unsigned int el = phi; el < n_subels(); el+=nPhi) {
		mvec vs = subelReaction(el, R);
		for(unsigned int i=0; i<vs.size(); i++)
			v[df_subindex(el,i)/nPhi] = vs[i];
	}
	return v;
}

void ReactiveUnitSet::_setDF(unsigned int DF, double v) {
	invert_index(DF, ixn_el, ixn_el_df);
	setSubelDF(ixn_el, ixn_el_df, v);
}
      
void ReactiveUnitSet::add_DF_group(unsigned int N) {
	for(unsigned int j=0; j<N; j++) {
		for(unsigned int i=0; i<nPhi; i++) {
			df_subel.push_back(n_subels()+i);
			df_subel_df.push_back(j);
		}
	}
	group_DF.push_back(N);
	group_start.push_back(group_start.back()+nPhi*N);
}

//----------------------------------------------

ReactiveSetCombiner::~ReactiveSetCombiner() {
	if(ownSets) {
		for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
			delete *it;
	}
	mySets.clear();
}

void ReactiveSetCombiner::addSet(ReactiveSet* R) {
	assert(R);
	assert(R->nPhi == nPhi);
	for(unsigned int i=0; i<R->nDF(); i++)
		df_set.push_back(mySets.size());
	mySets.push_back(R);
	set_cum_df.push_back(set_cum_df.back()+R->nDF());
}

void ReactiveSetCombiner::setInteractionDF(unsigned int DF, double v) {

	//std::cout << "RSC Setting Interaction DF " << DF << " of " << nDF() << " to " << v << std::endl;

	assert(DF < finalState.size() && ixn_df < finalState.size());
	finalState[ixn_df] = 0;
	finalState[DF] = v;
	
	ixn_set = df_set[ixn_df];
	mySets[ixn_set]->setInteractionDF(ixn_df-set_cum_df[ixn_set],0);
	
	ixn_set = df_set[DF];
	mySets[ixn_set]->setInteractionDF(DF-set_cum_df[ixn_set],v);
	
	ixn_df = DF;
}

void ReactiveSetCombiner::_setDF(unsigned int DF, double v) {
	assert(DF<nDF());
	unsigned int is = df_set[DF];
	mySets[is]->setDF(DF-set_cum_df[is],v);
}

void ReactiveSetCombiner::_setDFv(const mvec& v) {
	for(unsigned int i=0; i<mySets.size(); i++)
		mySets[i]->setFinalState(v.subvec(set_cum_df[i], set_cum_df[i+1]));
}

void ReactiveSetCombiner::startInteractionScan() {
	ReactiveSet::startInteractionScan();
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
		(*it)->startInteractionScan();
}

mvec ReactiveSetCombiner::getReactionTo(ReactiveSet* R, unsigned int phi) {
	mvec v(nDF()/nPhi);
	for(unsigned int i=0; i<mySets.size(); i++) {
		mvec vi = mySets[i]->getReactionTo(R,phi);
		v.load_subvec(vi, set_cum_df[i]/nPhi);
	}
	return v;
}

void ReactiveSetCombiner::prepareIncident() {
	if(incidentState.size() != nDF())
		incidentState = mvec(nDF());
	for(unsigned int i=0; i<mySets.size(); i++) {
		mySets[i]->prepareIncident();
		incidentState.load_subvec(mySets[i]->incidentState, set_cum_df[i]);
	}
}
