#include "ReactiveSet.hh"
#include <iostream>

void ReactiveSet::setInteractionDF(unsigned int DF, double v) {
	assert(DF < finalState.size() && ixn_df < finalState.size());
	finalState[ixn_df] = 0;
	finalState[DF] = v;
	ixn_df = DF;
}

//-----------------------------------------

/*
/// sub-element reaction to RS via protocol
	virtual mvec subelReaction() = 0;
	
	// interaction calculator cache variables
	mvec ic_v;			//< cached interaction matrix between sub-elements
	unsigned int ic_i;	//< i element interaction
	unsigned int ic_di;	//< i element DF
	unsigned int ic_DF;	//< reacting to overall DF
*/

mdouble ReactiveUnitSet::nextInteractionTerm(unsigned int& i, unsigned int& j) {
	// get previous cached term
	mdouble v = -ic_v[ic_di];	// NOTE - sign; TODO put this somewhere better
	i = df_subindex(ic_i, ic_di);
	j = ixn_df;
	
	// self-interaction TODO put this somewhere better
	if(i==j) v = 1.;
	
	// step to next term
	ic_di = (ic_di+1)%ic_v.size();
	if(!ic_di) {
		setInteractionDF((ixn_df+1)%nDF());
		if(!ixn_df) {
			ic_i = (ic_i+nPhi)%n_subels();
			if(ic_i < nPhi) ic_i = (ic_i+1)%nPhi;
		}
		ic_v = subelReaction();
	}
	
	return v;
}

void ReactiveUnitSet::setInteractionDF(unsigned int DF, double v) {
	invert_index(ixn_df, ixn_el, ixn_el_df);
	setSubelDF(ixn_el,ixn_el_df,0.);
	ReactiveSet::setInteractionDF(DF,v);
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

void ReactiveUnitSet::setFinalState(const mvec& v) {
	ReactiveSet::setFinalState(v);
	unsigned int el,df;
	for(unsigned int DF=0; DF<nDF(); DF++) {
		invert_index(DF,el,df);
		setSubelDF(el,df,v[DF]);
	}
}
