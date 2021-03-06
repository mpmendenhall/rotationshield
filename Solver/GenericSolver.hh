/* 
 * GenericSolver.hh, part of the RotationShield program
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

#ifndef GENERICSOLVER_HH
/// Make sure this header is only loaded once
#define GENERICSOLVER_HH

#include "Typedefs.hh"
#include "InteractionSolver.hh"
#include "gsl/gsl_matrix.h"
#ifdef WITH_LAPACKE
#include "LAPACKE_Matrix.hh"
#endif

/// Green's Function solver for systems of linear interactions

class GenericSolver: public InteractionSolver {
public:
    /// Constructor
    GenericSolver();
    /// Destructor
    virtual ~GenericSolver();
    /// Solve for the Greene's Function of a ReactiveSet system
    virtual void solve(ReactiveSet& R);
    /// Apply solution to ReactiveSet system initial state
    virtual void calculateResult(ReactiveSet& R);

#ifdef WITH_LAPACKE
    /// get un-normalized vector of singular values
    virtual mvec get_singular_values() const { if(my_SVD) return my_SVD->singular_values().getData(); return mvec(); }
#endif
        
protected:
    
    /// Assembles the interaction matrix
    void buildInteractionMatrix(ReactiveSet& R);

#ifdef WITH_LAPACKE
    mmat the_ixn;                               ///< The interaction matrix R between degrees of freedom
    LAPACKE_Matrix_SVD<double,double>*  my_SVD; ///< SVD of (1-R)
#else
    gsl_matrix* the_GF; ///< The inverted interaction matrix, i.e. the Green's Function for the system
#endif

};

#endif
