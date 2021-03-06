/* 
 * ReactiveElement.hh, part of the RotationShield program
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

#ifndef REACTIVEELEMENT_HH
/// Make sure this header is only loaded once
#define REACTIVEELEMENT_HH

#include "FieldSource.hh"
#include "Color.hh"
#include <cmath>
#include "VarMat.hh"

/// FieldSource that responds linearly to incident fields in order to satisfy boundary conditions
class ReactiveElement: public FieldSource {
public:
    
    /// constructor
    ReactiveElement() { }

    /// number of degrees of freedom
    virtual unsigned int nDF() const = 0;
    
    /// destructor
    virtual ~ReactiveElement() {}
        
    /// get fields at a point due to state components
    virtual mmat fieldAtComponents(vec3 p0) const = 0;
    /// resulting state in response to incident field f
    virtual mvec responseToFieldSource(const FieldSource* f) const = 0;
    /// interaction matrix with another ReactiveElement's state
    virtual mmat interactionWithComponents(const ReactiveElement* e) const = 0;
    /// interaction matrix with another ReactiveElement in its current state
    virtual mvec interactionWith(const ReactiveElement* e) const = 0;
    
    /// set state in reaction to applied field
    void reactTo(const FieldSource* f) { setState(responseToFieldSource(f)); }
    
    /// get field produced at a point
    vec3 fieldAt(const vec3& v) const { 
        mvec f =  fieldAtComponents(v)*state;
        return vec3(f[0],f[1],f[2]);
    }

    /// get state vector
    const mvec& getState() const { return state; }
    
    /// set state vector
    virtual void setState(mvec v) { assert(v.size() == nDF()); state = v; }
    /// set element of state vector
    virtual void setState(unsigned int df, double v) {
        if(state.size() != nDF()) state = mvec(nDF());
        assert(df<nDF());
        state[df] = v;
    }
        
protected:
    
    mvec state;         ///< vector describing state of element (e.g. surface current for a PlaneSource)
};

#endif
