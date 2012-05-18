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
	virtual void visualize(bool top = true, mdouble logmax = 3.0) const {
		
		if(top) { vsr::startRecording(true); vsr::clearWindow(); }
		
		float smag = state.mag2();
		if(smag)
			smag = max(0,min(1.0,0.1*(log(smag)+10-logmax)));
		
		vec3 svec = vec3();
		for(unsigned int i=0; i<min(3,state.size()); i++)
			svec[i] = state[i];
		
		vec3 hsv = vec3( atan2(svec[0],svec[1]), 1.0, 1.0 );
		vec3 rgb = hsv2rgb(hsv);
		vsr::setColor(rgb[0], rgb[1] , rgb[2], 0.5*smag);
		p.visualize(false);
		vsr::setColor(1, .7, .7, 1.0);
		p.visualizeCoordinates(0.2);
		vsr::setColor(.7, .7, 1, 1.0);
		p.visualizeCoordinates(-0.2);
		vsr::setColor(1,1,1,1);
		
		p.visualizeVector(svec);
		if(top) vsr::stopRecording();
	}
	
	Plane p;		//< The plane in which the element resides
};


#endif
