/// \file "geometry.hh" \brief Classes for describing shield geometry
#ifndef SHIELDBUILDER_HH
/// Make sure this file is only loaded once 
#define SHIELDBUILDER_HH 1

#include "MiscUtils.hh"
#include "RefCounter.hh"
#include "Geometry.hh"
#include "PlanarElement.hh"
#include <iostream>
#include <vector>
#include <cassert>
#include <math.h>

/// segment of a shield, generates a ring of PlanarElements
class ShieldSegment: public RefCounter {
public:
	ShieldSegment(unsigned int n, PlanarElement* b): nTheta(n), base(b) { base->retain(); }
	~ShieldSegment() { base->release(); }
	PlanarElement* genElement(unsigned int n) { return base->replicateRotated(2.0*n*M_PI/mdouble(nTheta)); }
	const unsigned int nTheta;
protected:
	
	PlanarElement* base;	//< PlanarElement that will be replicated into a ring
};


/// Rough estimate of field magnitude (in 2D plane) due to wires perpendicular to shield, used to optimize shield gridding for cos theta coil endcaps
class FieldEstimator2D {
public:
	/// Constructor
	FieldEstimator2D(): sources(std::vector<vec2>()), currents(std::vector<mdouble>()) {}
	/// Destructor
	virtual ~FieldEstimator2D() {}
	/// Estimated field at a point
	virtual vec2 estimateAt(const vec2& v) const;
	/// estimate rate of change in given direction
	virtual vec2 derivAt(const vec2& v, const vec2& dx) const;
	/// Add a "line source" perpendicular to plane of interest
	void addsource(const vec2& v, mdouble j) { sources.push_back(v); currents.push_back(j); }
protected:
	std::vector<vec2> sources;
	std::vector<mdouble> currents;
};

/// Field estimator based off 3D field source
class FieldEstimator2Dfrom3D: public FieldEstimator2D {
public:
	/// Constructor
	FieldEstimator2Dfrom3D(FieldSource* FS): myFS(FS) { assert(myFS); myFS->retain(); }
	/// Destructor
	virtual ~FieldEstimator2Dfrom3D() { myFS->release(); }
	/// Estimated field at a point
	virtual vec2 estimateAt(const vec2& v) const;
protected:
	FieldSource* myFS;
};


/// Base class for describing cylindrically-symmetric gridded surfaces (i.e. the magnetic shield)
class ShieldBuilder: public RefCounter {
public:
	/// Constructor
	ShieldBuilder(unsigned int nth): nTheta(nth), segments(std::vector<ShieldSegment*>()) {}
	/// Destructor
	~ShieldBuilder() { for(unsigned int i=0; i<nSegments(); i++) segments[i]->release(); }
	
	/// number of divisions along z axis
	unsigned int nSegments() const { return segments.size(); }
	/// access the Zth segment
	ShieldSegment* operator[](unsigned int z) const { assert(z<nSegments()); return segments[z]; }
	
	/// generate reactive surface element
	ReactiveElement* genElement(unsigned int z, unsigned int th) const { assert(z<nSegments()); return segments[z]->genElement(th); }
	
	/// add another ShieldBuilder's segments
	void append(const ShieldBuilder* G) { assert(G->nTheta == nTheta); for(unsigned int i=0; i<G->nSegments(); i++) append((*G)[i]); }
	/// add a ring of elements
	void append(ShieldSegment* s) { s->retain(); segments.push_back(s); }
	
	const unsigned int nTheta; //< number of divisions radially
	
	void OptCone(unsigned int nZ0, unsigned int nZ1, vec2 s, vec2 e, PlanarElement* base, FieldEstimator2D* fes = NULL);
	void makeOptCyl(unsigned int nZ0, unsigned int nZ1, mdouble r, mdouble z0, mdouble z1, PlanarElement* base, FieldEstimator2D* fe = NULL) { OptCone(nZ0,nZ1,vec2(z0,r),vec2(z1,r),base,fe); }
	
protected:
	std::vector<ShieldSegment*> segments;
};

#endif
