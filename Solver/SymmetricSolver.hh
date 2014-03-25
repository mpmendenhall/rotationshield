/// \file "SymmetricSolver.hh" \brief Contains the class for solving symmetric interacting systems
#ifndef SYMMETRICSOLVER_HH
/// Make sure this header is only loaded once
#define SYMMETRICSOLVER_HH 1

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "BlockCMat.hh"
#include "BinaryOutputObject.hh"
#include <string>
#include <cassert>

/// Green's Function solver for systems of linear interactions with a periodic symmetry between interaction terms (represented by ReactiveSets with nPhi > 1)

class SymmetricSolver: public InteractionSolver, public BinaryOutputObject {
public:
	/// Constructor
	SymmetricSolver(): InteractionSolver(true), singular_epsilon(1e-4), the_GF(NULL) {}
	/// Destructor
	virtual ~SymmetricSolver() { if(the_GF) delete the_GF; }
	
	/// Solve for the Greene's Function of a ReactiveSet system
	virtual void solve(ReactiveSet& R);
	/// Apply solution to ReactiveSet system initial state
	virtual void calculateResult(ReactiveSet& R);
	/// Apply self-interaction to final state
	virtual void selfInteract(ReactiveSet& R);
	
	/// Dump binary data to file
	void writeToFile(std::ostream& o) const;
	/// Read binary data from file
	static SymmetricSolver* readFromFile(std::istream& s);
	/// Read solution from file if available; otherwise, solve; save result to same file
	static SymmetricSolver* cachedSolve(ReactiveSet& R, const std::string& fname);

#ifdef WITH_LAPACKE
	/// set singular values threshold
	void set_singular_epsilon(double e) { singular_epsilon = e; }
#else
	/// set singular values threshold
	void set_singular_epsilon(double) { assert(false); }
#endif
	/// show singular values
	void print_singular_values() const;
	
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	/// perform in-place block circulant multiplication
	static void circulantMul(const BlockCMat& M, mvec& v, unsigned int nPhi);
	/// check inversion accuracy
	static double checkInversion(const BlockCMat& M, const BlockCMat& MI, unsigned int nPhi);
	
	BlockCMat the_ixn;			//< The interaction matrix R between degrees of freedom
	double singular_epsilon;	//< singular value threshold
	BlockCMat_SVD* the_GF;		//< the Green's Function for the system, (I-R)^-1
};


#endif
