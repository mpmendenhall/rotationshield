/* 
 * ReactiveSet.cpp, part of the RotationShield program
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

#include "ReactiveSet.hh"
#include "ProgressBar.hh"
#include <iostream>

void ReactiveSet::setDF(unsigned int DF, double v) {
	if(finalState.size() != nDF())
		finalState = mvec(nDF());
	finalState[DF] = v;
	_setDF(DF,v);
}

void ReactiveSet::setInteractionDF(unsigned int DF, double v) {
	if(DF < nDF()) {
		if(ixn_df < nDF()) {
			finalState[ixn_df] = 0;
			_setDF(ixn_df,0);
		}
		finalState[DF] = v;
		_setDF(DF,v);
	}
	ixn_df = DF;
}

void ReactiveSet::_setDFv(const mvec& v) {
	assert(v.size() == nDF());
	for(unsigned int DF = 0; DF < nDF(); DF++)
		_setDF(DF,v[DF]);
}

mvec ReactiveSet::getFullReactionTo(ReactiveSet* R) {
	printf("Calculating reaction of %i = %i x %i elements...\n", nDF(), nDF()/nPhi, nPhi);
	ProgressBar pb = ProgressBar(nPhi, 1, true);
	
	mvec v(nDF());
	for(unsigned int p=0; p<nPhi; p++) {
		mvec vp = getReactionTo(R,p);
		for(unsigned int i=0; i<vp.size(); i++)
			v[i*nPhi+p] = vp[i];
		pb.update(p);
	}
	return v;
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
	ReactiveSet::setInteractionDF(DF,v);
	if(DF >= nDF()) {
		for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
			(*it)->setInteractionDF((*it)->nDF(),0);
	} else {
		ixn_set = df_set[DF];
	}
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

