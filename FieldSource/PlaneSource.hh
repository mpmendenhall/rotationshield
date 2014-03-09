#ifndef PLANESOURCE_HH
#define PLANESOURCE_HH 1

#include "PlanarElement.hh"
#include "Matrix.hh"

/// Magnetic field source due to a plane of current (e.g. a bound surface current on a magnetic shield)
class PlaneSource: public PlanarElement {
public:
	/// Constructor
	PlaneSource(Plane cp, mdouble mu): PlanarElement(cp), murel(mu) { setRmat(); }
	
	/// number of degrees of freedom
	virtual unsigned int nDF() const { return 2; }
	
	/// state in response to applied field
	virtual mvec responseToFieldSource(const FieldSource* f) const { return rmat * f->fieldOverPlane(p); }
	
	/// interaction matrix with another ReactiveElement
	virtual mmat interactionWithComponents(const ReactiveElement* e) const { return rmat * e->fieldAtComponents(p.o); }
	/// interaction matrix with another ReactiveElement
	virtual mvec interactionWith(const ReactiveElement* e) const { return rmat * e->fieldAt(p.o); }

	/// replicate around a new plane
	virtual PlanarElement* replicate(Plane pl) const { return new PlaneSource(pl,murel); }
	
	/// generate initial reference element
	virtual PlanarElement* reference(annulusSpec a) const { return new PlaneSource(Plane(a),murel); }
	
	/// replicate reference element to other angle
	virtual PlanarElement* replicateRotated(mdouble th) const { return new PlaneSource(p.zrotated(th),murel); }

	/// field components due to each DF at given point
	virtual mmat fieldAtComponents(vec3 p0) const;
			
private:
	mdouble murel;	//< relative permeability
	mmat rmat;		//< response matrix to applied field
	
	/// generate correct response matrix to applied fields
	void setRmat() {
		Matrix<2,3,mdouble> M = Matrix<2,3,mdouble>();
		M(0,1) = M(1,0) = 2.0*(1.0-murel)/(1.0+murel);
		rmat = mmat(M*p.projectionMatrix());		
	}
	
};

#endif
