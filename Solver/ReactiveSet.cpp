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

mdouble ReactiveUnitSet::nextInteractionTerm(unsigned int& i, unsigned int& j) {
	// get previous cached term
	mdouble v = ic_v[ic_di];
	i = df_subindex(ic_i, ic_di);
	j = ixn_df;
	
	// self-interaction
	if(i==j) v = 0;
	
	// step to next term
	ic_di = (ic_di+1)%ic_v.size();
	if(!ic_di) {
		setInteractionDF((ixn_df+1)%nDF());
		invert_index(ixn_df, ixn_el, ixn_el_df);
		if(!ixn_df) {
			ic_i = (ic_i+nPhi)%n_subels();
			if(ic_i < nPhi) ic_i = (ic_i+1)%nPhi;
		}
		ic_v = subelReaction();
	}
	
	return v;
}

void ReactiveUnitSet::_setDF(unsigned int DF, double v) {
	unsigned int el, eldf;
	invert_index(DF, el, eldf);
	setSubelDF(el, eldf, v);
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

bool ReactiveSetCombiner::set_protocol(void* ip) {
	ReactiveUnitSet::set_protocol(ip);
	bool b = true;
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
		b &= (*it)->set_protocol(ip);
	return b;
}

void ReactiveSetCombiner::setSubelDF(unsigned int el, unsigned int df, mdouble v) {
	assert(false); //TODO
	//mySets[i]->
}

void ReactiveSetCombiner::setInteractionDF(unsigned int DF, double v=1.0) {
	
}

void ReactiveSetCombiner::setFinalState(const mvec& v) {
	assert(v.size()==nDF());
	finalState = v;
	for(unsigned int i=0; i<mySets.size(); i++)
		mySets[i]->setFinalState(finalState.subvec(group_start[i],group_start[i+1]));
}

void ReactiveSetCombiner::setZeroState() {
	ReactiveUnitSet::setZeroState();
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
		(*it)->setZeroState();
}

*/

