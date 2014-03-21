#ifndef GENERICSOLVER_HH
#define GENERICSOLVER_HH 1

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "gsl/gsl_matrix.h"

/// Green's Function solver for systems of linear interactions

class GenericSolver: public InteractionSolver {
public:
	/// Constructor
	GenericSolver(): InteractionSolver(true), the_GF(NULL) {}
	/// Destructor
	virtual ~GenericSolver();
	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R);
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R);
		
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	gsl_matrix* the_GF; //< The inverted interaction matrix, i.e. the Green's Function for the system
};

#endif
