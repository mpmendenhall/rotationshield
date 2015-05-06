/* 
 * MultiQuilibrator.hh, part of the RotationShield program
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

 /// \file "MultiQuilibrator.hh" \brief Bring several individual solved systems into mutual equilibrium
#ifndef MULTIQUILIBRATOR_HH
/// Make sure this header is only loaded once
#define MULTIQUILIBRATOR_HH

#include "InteractionSolver.hh"
#include <vector>

struct solvedSet {
    ReactiveSet* RS;            ///< the interacting system
    InteractionSolver* IS;      ///< its solver
    mvec v_inc;                 ///< initial incident response
    mvec v_ext;                 ///< sequence of perturbation responses
};

/// Iteratively equilibrates multiple ReactiveSets
class MultiQuilibrator {
public:
    /// constructor
    MultiQuilibrator() {}
    /// destructor
    virtual ~MultiQuilibrator() {}
    
    /// add a subset
    void addSet(ReactiveSet* RS, InteractionSolver* IS);
    
    /// perform one equilibration step
    double step();
    /// continue stepping until relative error limit reached
    unsigned int equilibrate(double rel_err = 1e-4);
    
protected:

    /// update one set
    double update_set(unsigned int i);
    
    std::vector<solvedSet> mySets;      ///< interacting systems
    
};

#endif
