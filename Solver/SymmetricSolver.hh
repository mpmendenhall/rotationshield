/// \file "SymmetricSolver.hh" \brief Contains the class for solving symmetric interacting systems
#ifndef SYMMETRICSOLVER_HH
/// Make sure this header is only loaded once
#define SYMMETRICSOLVER_HH 1

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "BlockCMat.hh"
#include <string>

/// Green's Function solver for systems of linear interactions with a periodic symmetry between interaction terms (represented by ReactiveSets with nPhi > 1)

class SymmetricSolver: public InteractionSolver {
public:
	/// Constructor
	SymmetricSolver(): InteractionSolver(true) {}

	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R);
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R);
	
	/// Write solution to file
	void writeToFile(std::ostream& o) const;
	/// Read solution from file
	void readFromFile(std::istream& s) const;
	/// Read solution from file if available; otherwise, solve; save result to same file
	void cachedSolve(ReactiveSet& R, const std::string& fname);
	
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	BlockCMat<mdouble> the_GF; //< The inverted interaction matrix, i.e. the Green's Function for the system
};


#endif
