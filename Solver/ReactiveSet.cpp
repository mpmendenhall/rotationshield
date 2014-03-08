#include "ReactiveSet.hh"
#include <iostream>

mdouble ReactiveUnitSet::nextInteractionTerm(unsigned int& i, unsigned int& j) {
	// get previous cached term
	mdouble v = -ic_m(ic_di, ic_dj);	// NOTE - sign; TODO put this somewhere better
	i = df_subindex(ic_i, ic_di);
	j = df_subindex(ic_j, ic_dj);
	
	// self-interaction TODO put this somewhere better
	//if(ic_i==ic_j) v = (ic_di==ic_dj);
	if(ic_i==ic_j && ic_di==ic_dj) v = 1;
	
	// step to next term
	ic_dj = (ic_dj+1)%ic_m.nCols();
	if(!ic_dj) {
		ic_di = (ic_di+1)%ic_m.nRows();
		if(!ic_di) {
		
			ic_j = (ic_j+1)%n_subels();
			if(!ic_j) {
				ic_i = (ic_i+nPhi)%n_subels();
				if(ic_i < nPhi) ic_i = (ic_i+1)%nPhi;
			}
			
			ic_m = interactionBetween(ic_i, ic_j);
		}
	}

	return v;
}
