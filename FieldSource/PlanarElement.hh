#ifndef PLANARELEMENT_HH
#define PLANARELEMENT_HH 1

#include "ReactiveElement.hh"

/// a ReactiveElement on a plane
class PlanarElement: public ReactiveElement {
public:
	/// constructor
	PlanarElement(Plane pl): ReactiveElement(), p(pl) { }
		
	/// replicate around a new plane
	virtual PlanarElement* replicate(Plane pl) const = 0;
	
	/// replicate with new annulus specification
	virtual PlanarElement* reference(annulusSpec a) const { return replicate(Plane(a)); }
	
	/// replicate rotated around z axis
	virtual PlanarElement* replicateRotated(mdouble th) const = 0;
	
	/// Visualize the element
	virtual void visualize(bool top = true, mdouble logmax = 3.0) const;
	
	Plane p;		//< The plane in which the element resides
};


#endif
