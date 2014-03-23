/// \file "SymmetricSolver.hh" \brief Contains the class for solving symmetric interacting systems
#ifndef SYMMETRICSOLVER_HH
/// Make sure this header is only loaded once
#define SYMMETRICSOLVER_HH 1

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "CMatrix.hh"
#include "VarMat.hh"
#include <string>

typedef  VarMat<CMatrix> BlockCMat;

/// Green's Function solver for systems of linear interactions with a periodic symmetry between interaction terms (represented by ReactiveSets with nPhi > 1)

class SymmetricSolver: public InteractionSolver {
public:
	/// Constructor
	SymmetricSolver(): InteractionSolver(true) {}

	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R);
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R);
	/// Apply self-interaction to final state
	virtual void selfInteract(ReactiveSet& R);
	
	/// Write solution to file
	void writeToFile(std::ostream& o) const;
	/// Read solution from file
	void readFromFile(std::istream& s) const;
	/// Read solution from file if available; otherwise, solve; save result to same file
	void cachedSolve(ReactiveSet& R, const std::string& fname);
	
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	/// perform in-place block circulant multiplication
	static void circulantMul(const BlockCMat& M, mvec& v, unsigned int nPhi);
	/// check inversion accuracy
	static double checkInversion(const BlockCMat& M, const BlockCMat& MI, unsigned int nPhi);
	/// construct block circulant identity
	static BlockCMat makeIdentity(unsigned int N, unsigned int nPhi);
	
	BlockCMat the_ixn;	//< The interaction matrix R between degrees of freedom
	BlockCMat the_GF;	//< the Green's Function for the system, (I-R)^-1
};


#endif
