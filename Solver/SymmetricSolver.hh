/// \file "SymmetricSolver.hh" \brief Contains the core class for doing shielding simulations
#ifndef SYMMETRICSOLVER_HH
/// Make sure this header is only loaded once
#define SYMMETRICSOLVER_HH 1

#include "MiscUtils.hh"
#include "InteractionSolver.hh"
#include "BlockCMat.hh"

/// The class where everything else is pulled together to solve the shield boundary-value problem.

/// Usage:
///		- Create a SurfacelCyl describing the geometry of the shield
///		- Create an instance of the class with SymmetricSolver(SurfacelCyl* g)
///		- Call solve() to calculate the shield's Green's Function (this takes time... we have to calculate and invert the interaction matrix)
///		- Create a FieldSource describing the incident field that you want to see the shield's response to
///		- Call calculate_incident() using the FieldSource (this takes time, depending on the field source)
///		- Call calculateResult()
///		- Use a FieldAnalyzer to record the resulting magnetic fields

class SymmetricSolver: public InteractionSolver {
public:
	/// Constructor; uses SurfacelCyl * g to describe the shield geometry
	SymmetricSolver();
	/// Destructor
	virtual ~SymmetricSolver() {}
	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R);
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R);
		
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	BlockCMat<mdouble> the_GF; //< The inverted interaction matrix, i.e. the Green's Function for the system
};


#endif
