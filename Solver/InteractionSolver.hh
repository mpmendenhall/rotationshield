/* 
 * InteractionSolver.hh, part of the RotationShield program
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

#ifndef INTERACTIONSOLVER_HH
#define INTERACTIONSOLVER_HH 1

#include "ReactiveSet.hh"
#include "Typedefs.hh"
#include <stdlib.h>
#include <math.h>

/// virtual base class for solving interacting systems
class InteractionSolver {
public:
	/// Constructor
	InteractionSolver(bool vbs = true): verbose(vbs) {}
	/// Destructor
	virtual ~InteractionSolver() {}

	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R) = 0;
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R) = 0;
	
	bool verbose;	//< whether to display solver progress
};

#endif
