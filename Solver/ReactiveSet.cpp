#include "ReactiveSet.hh"
#include <iostream>

void ReactiveSet::setInteractionDF(unsigned int DF, double v) {
	assert(DF < finalState.size() && ixn_df < finalState.size());
	finalState[ixn_df] = 0;
	finalState[DF] = v;
	_setDF(ixn_df,0);
	_setDF(DF,v);
	ixn_df = DF;
}

void ReactiveSet::_setDF(const mvec& v) {
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

/*

void ReactiveSetCombiner::addSet(ReactiveSet* R) {
	assert(R->nPhi == nPhi);
	mySets.push_back(R);
	add_DF_group(R->nDF()/nPhi);
}

void ReactiveSetCombiner::setSubelDF(unsigned int el, unsigned int df, mdouble v) {
	mySets[el]->setDF(df,v);
}

void ReactiveSetCombiner::_setDF(const mvec& v) {
	for(unsigned int i=0; i<mySets.size(); i++)
		mySets[i]->setDF(v.subvec(group_start[i], group_start[i+1]));
}

void ReactiveSetCombiner::startInteractionScan() {
	ReactiveUnitSet::startInteractionScan();
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
		(*it)->startInteractionScan();
}

*/