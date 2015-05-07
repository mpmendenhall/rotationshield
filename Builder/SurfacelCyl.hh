/* 
 * SurfacelCyl.hh, part of the RotationShield program
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

/// \file SurfacelCyl.hh \brief Cylindrical arrangements of surface elements (replacement for legacy SurfacelCyl)

#ifndef SurfacelCyl_HH
/// Make sure this header is only loaded once
#define SurfacelCyl_HH

#include "ReactiveSet.hh"
#include "ReactiveElement.hh"
#include "FieldSource.hh"
#include "FieldEstimator2D.hh"
#include "PlanarElement.hh"
#include "MagRS.hh"

/// surface element interaction protocol
class Surfacel_Protocol {
public:
    ReactiveElement* e;
    static Surfacel_Protocol* SP;
};

/// collection of individual surface elements
class SurfacelSet: public ReactiveUnitSet, public FieldSource {
public:
    /// constructor
    SurfacelSet(unsigned int nph): ReactiveUnitSet(nph), FieldSource(), verbose(true) { }
    /// destructor
    virtual ~SurfacelSet();
    
    /// get field produced at a point
    virtual vec3 fieldAt(const vec3& v) const;
    
    /// add a surface element
    void addSurfacel(ReactiveElement* e);
    
    /// visualize surface elements
    virtual void _visualize() const;
    
    /// calculate response to incident field
    virtual void calculateIncident(const FieldSource& f);
    
    /// respond to interaction protocol; respond if protocol identified
    virtual bool queryInteraction(void* ip);
    
    bool verbose;    ///< whether to display calculation progress
    
protected:
    vector<ReactiveElement*> surfacels;    ///< the surface elements
    
    //======================================
    /// set state for i^th sub-element
    virtual void setSubelDF(unsigned int el, unsigned int df, double v) { surfacels[el]->setState(df,v); }
    /// sub-element reaction to RS via protocol
    virtual mvec subelReaction(unsigned int el, ReactiveSet* R);
    //======================================
};



//--------------------------------------------------
//--------------------------------------------------
//--------------------------------------------------


/// segment of a shield, generates a ring of PlanarElements
class ShieldSegment: public RefCounter {
public:
    ShieldSegment(unsigned int n, PlanarElement* b): nTheta(n), base(b) { base->retain(); }
    ~ShieldSegment() { base->release(); }
    PlanarElement* genElement(unsigned int n) { return base->replicateRotated(2.0*n*M_PI/double(nTheta)); }
    const unsigned int nTheta;
protected:
    
    PlanarElement* base;    ///< PlanarElement that will be replicated into a ring
};

/// Base class for describing cylindrically-symmetric gridded surfaces (i.e. the magnetic shield)
class SurfacelCyl: public SurfacelSet {
public:
    /// Constructor
    SurfacelCyl(unsigned int nth): SurfacelSet(nth) { }
    /// Destructor
    ~SurfacelCyl();
            
    /// add another SurfacelCyl's segments
    void append(const SurfacelCyl* G) { assert(G->nPhi == nPhi); for(unsigned int i=0; i<G->segments.size(); i++) append(G->segments[i]); }
    
    /// add a ring of elements
    void append(ShieldSegment* s);
    
    /// make a general conical shield
    void OptCone(unsigned int nZ0, unsigned int nZ1, vec2 s, vec2 e, PlanarElement* base, FieldEstimator2D* fes = NULL);
    
    /// make a cylindrical shield
    void makeOptCyl(unsigned int nZ0, unsigned int nZ1, double r, double z0, double z1, PlanarElement* base, FieldEstimator2D* fe = NULL) { OptCone(nZ0,nZ1,vec2(z0,r),vec2(z1,r),base,fe); }
    
protected:
    vector<ShieldSegment*> segments;    ///< element generators for each ring in shield
};



#endif
