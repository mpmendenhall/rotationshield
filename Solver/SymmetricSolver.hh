/* 
 * SymmetricSolver.hh, part of the RotationShield program
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

/// \file "SymmetricSolver.hh" \brief Contains the class for solving symmetric interacting systems
#ifndef SYMMETRICSOLVER_HH
/// Make sure this header is only loaded once
#define SYMMETRICSOLVER_HH

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
	SymmetricSolver(): InteractionSolver(true), the_GF(NULL) {}
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
	/// get un-normalized vector of singular values
	virtual mvec get_singular_values() const;
	/// get right singular vector
	virtual mvec get_singular_vector(unsigned int i) const;
#endif
	
protected:
	
	/// Assembles the interaction matrix
	void buildInteractionMatrix(ReactiveSet& R);
	
	/// perform in-place block circulant multiplication
	static void circulantMul(const BlockCMat& M, mvec& v, unsigned int nPhi);
	/// check inversion accuracy
	static double checkInversion(const BlockCMat& M, const BlockCMat& MI, unsigned int nPhi);
	
	BlockCMat the_ixn;			///< The interaction matrix R between degrees of freedom
	BlockCMat_SVD* the_GF;		///< the Green's Function for the system, (I-R)^-1
};


#endif
