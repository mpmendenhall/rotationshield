#ifndef GENERICSOLVER_HH
#define GENERICSOLVER_HH 1

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "gsl/gsl_matrix.h"

/// The class where everything else is pulled together to solve the shield boundary-value problem.

/// Usage:
///		- Create a SurfacelCyl describing the geometry of the shield
///		- Create an instance of the class with GenericSolver(SurfacelCyl* g)
///		- Call solve() to calculate the shield's Green's Function (this takes time... we have to calculate and invert the interaction matrix)
///		- Create a FieldSource describing the incident field that you want to see the shield's response to
///		- Call calculate_incident() using the FieldSource (this takes time, depending on the field source)
///		- Call calculateResult()
///		- Use a FieldAnalyzer to record the resulting magnetic fields

class GenericSolver: public InteractionSolver {
public:
	/// Constructor; uses SurfacelCyl * g to describe the shield geometry
	GenericSolver();
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
