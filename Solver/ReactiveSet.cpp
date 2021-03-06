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

unsigned int ReactiveSet::n_reactive_sets = 0;

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
    bool rsym = R->isRotSym();
    if(rsym) printf("\tUsing rotational symmetry shortcut!\n");
    ProgressBar pb = ProgressBar(nPhi, 1, true);
    
    mvec v(nDF());
    for(unsigned int p=0; p<nPhi; p++) {
        if(p && rsym) {
            for(unsigned int i=0; i<nDF()/nPhi; i++)
                v[i*nPhi+p] = v[i*nPhi];
        } else {
            mvec vp = getReactionTo(R,p);
            for(unsigned int i=0; i<vp.size(); i++)
                v[i*nPhi+p] = vp[i];
        }
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

//-----------------------------------------



// Sub-DF ordering example:
//----------------------
//    6    3    2    1    nPhi
// ___ ___ ___ ___
//    a1    a    a  _a_
//    b1    c  _d_ _b_
//    c1 _e_    b  _c_
//    d1    b  _e_ _d_
//    e1    d    c  _e_
// _f1 _f_ _f_ _f_
// ___ ___ ___ ___
//    a2    a    a  _a_
//    b2    c  _d_ _b_
//    c2 _e_    b  _c_
//    d2    b  _e_ _d_
//    e2    d    c  _e_
// _f2 _f_ _f_ _f_
//
// getReactionTo(R, phi):
// 6: [a1,a2], ... [f1,f2]
// 3: [a1,b1,a2,b2], ... [e1,f1,e2,f2]


ReactiveSetCombiner::~ReactiveSetCombiner() {
    if(ownSets) {
        for(vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
            delete *it;
    }
    mySets.clear();
}

bool ReactiveSetCombiner::queryInteraction(void* ip) {
    if(interactionMode()) {
        assert(ixn_set < mySets.size());
        return mySets[ixn_set]->queryInteraction(ip);
    }
    bool did_interact  = false;
    for(vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
        did_interact |= (*it)->queryInteraction(ip);
    return did_interact;
}

void ReactiveSetCombiner::addSet(ReactiveSet* R) {
    assert(R);
    assert(!(R->nPhi%nPhi));
    for(unsigned int i=0; i<R->nDF(); i++) {
        df_set.push_back(mySets.size());
        df_set_df.push_back( R->nPhi*(i/R->nPhi) + (i%nPhi)*(R->nPhi/nPhi) + (i%R->nPhi)/nPhi );
    }
    mySets.push_back(R);
}

void ReactiveSetCombiner::setInteractionDF(unsigned int DF, double v) {
    if(DF >= nDF()) {
        for(vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
            (*it)->setInteractionDF((*it)->nDF(),0);
    } else {
        if(ixn_df < nDF()) {
            finalState[ixn_df] = 0;
            ixn_set = df_set[ixn_df];
            mySets[ixn_set]->setInteractionDF(df_set_df[ixn_df],0);
        }
        
        finalState[DF] = v;
        ixn_set = df_set[DF];
        mySets[ixn_set]->setInteractionDF(df_set_df[DF],v);
    }
    
    ixn_df = DF;
}

void ReactiveSetCombiner::_setDF(unsigned int DF, double v) {
    assert(DF<nDF());
    mySets[df_set[DF]]->setDF(df_set_df[DF],v);
}

void ReactiveSetCombiner::startInteractionScan() {
    ReactiveSet::startInteractionScan();
    for(vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++)
        (*it)->startInteractionScan();
}

mvec ReactiveSetCombiner::getReactionTo(ReactiveSet* R, unsigned int phi) {
    mvec v(nDF()/nPhi);
    unsigned int cum_df = 0;
    for(unsigned int i=0; i<mySets.size(); i++) {
        for(unsigned int j=0; j<(mySets[i]->nPhi/nPhi); j++) {
            mvec vi = mySets[i]->getReactionTo(R, phi*(mySets[i]->nPhi/nPhi) + j);
            for(unsigned int k=0; k<vi.size(); k++)
                v[cum_df + k*(mySets[i]->nPhi/nPhi) + j] = vi[k];
        }
        cum_df += mySets[i]->nDF()/nPhi;
    }
    return v;
}

void ReactiveSetCombiner::load_component_DF() {
    finalState = mvec(nDF());
    for(unsigned int DF=0; DF<nDF(); DF++)
        finalState[DF] = mySets[df_set[DF]]->finalState[df_set_df[DF]];
}


