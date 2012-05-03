#ifndef GENERICSOLVER_HH
#define GENERICSOLVER_HH 1

#include "MiscUtils.hh"
#include "InteractionSolver.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"

/// The class where everything else is pulled together to solve the shield boundary-value problem.

/// Usage:
///		- Create a ShieldBuilder describing the geometry of the shield
///		- Create an instance of the class with GenericSolver(ShieldBuilder* g)
///		- Call solve() to calculate the shield's Green's Function (this takes time... we have to calculate and invert the interaction matrix)
///		- Create a FieldSource describing the incident field that you want to see the shield's response to
///		- Call calculate_incident() using the FieldSource (this takes time, depending on the field source)
///		- Call calculateResult()
///		- Use a FieldAnalyzer to record the resulting magnetic fields

class GenericSolver: public InteractionSolver {
public:
	/// Constructor; uses ShieldBuilder * g to describe the shield geometry
	GenericSolver(): InteractionSolver() {}
	/// Destructor
	virtual ~GenericSolver() {}
	/// Solves the boundary value problem for the shield's geometry, producing #the_GF
	virtual void solve(bool checkinversion = false);
	/// Calculate the shield's #finalState in response to the #incidentState
	virtual void calculateResult();
	
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix();
	
	CLHEP::HepMatrix the_GF; //< The inverted interaction matrix, i.e. the Green's Function for the system
};

#endif
