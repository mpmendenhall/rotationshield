/// \file SurfaceSource.hh \brief Base class for extended 2D surface magnetic field sources

#ifndef SURFACESOURCE_HH
/// Makes sure to only load this file once
#define SURFACESOURCE_HH 1

#include "FieldSource.hh"
#include "Integrator.hh"

class SurfaceSource: public FieldSource {
public:
	/// destructor
	virtual ~SurfaceSource() {}
	
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, mdouble x, mdouble y) const = 0x0;
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const;
	
protected:
	mutable Integrator2D myIntegrator;	//< surface field integrator
};


#endif
