/* 
 * MultiQuilibrator.cpp, part of the RotationShield program
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

#include "MultiQuilibrator.hh"
#include <cassert>
#include <cmath>

void MultiQuilibrator::addSet(ReactiveSet* RS, InteractionSolver* IS) {
	assert(RS && IS);
	solvedSet s;
	s.RS = RS;
	s.IS = IS;
	mySets.push_back(s);
}

double MultiQuilibrator::update_set(unsigned int i) {
	assert(i<mySets.size());
	mySets[i].RS->incidentState = mySets[i].v_inc + mySets[i].v_ext;
	mvec prev = mySets[i].RS->finalState;
	mySets[i].IS->calculateResult(*mySets[i].RS);
	double di = (prev-mySets[i].RS->finalState).mag2();
	double mi = prev.mag2();
	std::cout << "Set " << i << " changed by " << sqrt(di) << " (" << 100*sqrt(di/mi) << " %)\n";
	return di/mi;
}

double MultiQuilibrator::step() {
	
	std::cout << "Equilibration step on " << mySets.size() << " systems...\n";
	double dv = 0;
	for(unsigned int i=0; i<mySets.size(); i++) {
	
		if(!mySets[i].v_inc.size()) mySets[i].v_inc = mySets[i].RS->incidentState;
		if(!mySets[i].v_ext.size()) mySets[i].v_ext = mvec(mySets[i].RS->nDF());
		
		// accumulate reaction to other systems
		mySets[i].v_ext.zero();
		for(unsigned int j=0; j<mySets.size(); j++) {
			if(i==j) continue;
			mySets[i].v_ext += mySets[i].RS->getFullReactionTo(mySets[j].RS);
		}
		
		dv += update_set(i);
	}
	
	std::cout << "\n** Net change\t" << 100*sqrt(dv) << "%\n\n";
	return sqrt(dv);
}

unsigned int MultiQuilibrator::equilibrate(double rel_err) {
	unsigned int nsteps = 0;
	while(step() > rel_err) { nsteps++; }
	return nsteps;
}
