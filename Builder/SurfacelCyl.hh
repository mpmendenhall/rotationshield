/// \file SurfacelCyl.hh \brief Cylindrical arrangements of surface elements (replacement for legacy SurfacelCyl)

#ifndef SurfacelCyl_HH
/// Makes sure to only load this file once
#define SurfacelCyl_HH 1

#include "ReactiveSet.hh"
#include "ReactiveElement.hh"
#include "FieldSource.hh"
#include "FieldEstimator2D.hh"
#include "PlanarElement.hh"

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
	SurfacelSet(unsigned int nph): ReactiveUnitSet(nph), FieldSource(), verbose(true) { set_protocol(Surfacel_Protocol::SP); }
	/// destructor
	virtual ~SurfacelSet();
	
	/// get field produced at a point
	virtual vec3 fieldAt(const vec3& v) const;
	
	/// add a surface element
	void addSurfacel(ReactiveElement* e);
	
	/// visualize surface elements
	virtual void visualize(bool top = true, mdouble scale = 1.0) const;
	
	/// calculate response to incident field
	virtual void calculateIncident(FieldSource* f);
	
	/// set state for i^th sub-element
	virtual void setState(unsigned int i, const mvec& v) { surfacels[i]->setState(v); }
	
	/// set interaction protocol to use; respond whether accepted
	virtual bool set_protocol(void* ip);
	/// respond to interaction protocol
	virtual void queryInteraction();
	
	bool verbose;	//< whether to display calculation progress
	
protected:
	std::vector<ReactiveElement*> surfacels;	//< the surface elements
	
	/// set state for i^th sub-element
	virtual void setSubelDF(unsigned int el, unsigned int df, mdouble v) { surfacels[el]->setState(df,v); }
	
	/// sub-element reaction to RS via protocol
	virtual mvec subelReaction();
};

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
	void makeOptCyl(unsigned int nZ0, unsigned int nZ1, mdouble r, mdouble z0, mdouble z1, PlanarElement* base, FieldEstimator2D* fe = NULL) { OptCone(nZ0,nZ1,vec2(z0,r),vec2(z1,r),base,fe); }
	
protected:
	std::vector<ShieldSegment*> segments;	//< element generators for each ring in shield
};



#endif
