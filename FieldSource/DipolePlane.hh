#ifndef DIPOLEPLANE_HH
#define DIPOLEPLANE_HH 1

#include "ReactiveElement.hh"
#include "PlaneSource.hh"

/// Magnetic field source due to a plane of dipoles (e.g. a thin shield)
class DipolePlane: public ReactiveElement {
public:
	/// Constructor
	DipolePlane(Plane cp, mdouble mu, mdouble thickness):  ReactiveElement(cp), d(thickness), murel(mu) { genSurfaces(); }
	
	/// Constructor, annular form
	DipolePlane(annulusSpec a, mdouble mu, mdouble thickness): ReactiveElement(a), d(thickness), murel(mu) { genSurfaces(); }
	
	/// Destructor
	virtual ~DipolePlane() { for(unsigned int i=0; i<6; i++) surfaces[i]->release(); }
	
	/// replicate around a new plane
	virtual ReactiveElement* replicate(Plane pl) const { return new DipolePlane(pl,murel,d); }
	
	/// generate initial reference element
	virtual ReactiveElement* reference(annulusSpec a) const { return new DipolePlane(a,murel,d); }
	
	/// replicate reference element to other angle
	virtual ReactiveElement* replicateRotated(mdouble th) const { return new DipolePlane(p.zrotated(th), murel, d); }
	
	/// Fields due to x- and z-directed components of the current
	mat3 fieldAtComponents(vec3 p0) const;
	
	/// nearfield response
	mat3 fieldAtComponentsNearfield(vec3 p0) const;
	
	/// Visualize the element
	virtual void visualize(bool top = true, mdouble logmax = 3.0) const {
		if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); }
		for(unsigned int i=0; i<drawSurfaces; i++)
			surfaces[i]->visualize(false);
		if(top) vsr::Visr::W->stopRecording();
	}
	
	/// point at which to measure incident field
	virtual vec3 measurePoint() const { return p.o + p.sn*0.290*d; } //surfaces[1]->p.o; }
	
	vec3 fieldAt(const vec3& v) const { vec3 f = vec3(); for(unsigned int i=0; i<nSurfaces; i++) f+=surfaces[i]->fieldAt(v); return f; }
	
	virtual void setState(vec3 v) {
		state = v;
		v = vec3(v[1],-v[0],0)/(2*d);
		surfaces[0]->setState(v);
		surfaces[1]->setState(v);
		surfaces[2]->setState(vec3(0,v[0],0));
		surfaces[3]->setState(vec3(0,v[0],0));
		surfaces[4]->setState(vec3(0,v[1],0));
		surfaces[5]->setState(vec3(0,v[1],0));
	}
	
	static bool doNearfield;
	static bool alwaysDipole;
	static unsigned int drawSurfaces;
	static unsigned int nSurfaces;
	
private:
	vec2 kfactor;
	mdouble d;
	mdouble murel;
	PlaneSource* surfaces[6];
	void genSurfaces();
};

#endif
