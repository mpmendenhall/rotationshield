#ifndef INTERACTIONSOLVER_HH
#define INTERACTIONSOLVER_HH 1

#include "ReactiveSet.hh"
#include "MiscUtils.hh"
#include <stdlib.h>
#include <math.h>

/// virtual base class for solving interacting systems
class InteractionSolver {
public:
	/// Constructor
	InteractionSolver() {}
	/// Destructor
	virtual ~InteractionSolver() {}

	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R) = 0;
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R) = 0;
	
	bool verbose;	//< whether to display solver progress
};

#endif
