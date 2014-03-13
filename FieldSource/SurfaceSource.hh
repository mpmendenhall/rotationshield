/// \file SurfaceSource.hh \brief Base class for extended 2D surface magnetic field sources

#ifndef SURFACESOURCE_HH
/// Makes sure to only load this file once
#define SURFACESOURCE_HH 1

#include "FieldSource.hh"
#include "Integrator.hh"

class SurfaceSource: public FieldSource {
public:
	/// constructor
	SurfaceSource(): FieldSource(), dflt_integrator_ndivs_x(1), dflt_integrator_ndivs_y(1) {}

	/// destructor
	virtual ~SurfaceSource() {}
	
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, const vec2& l) const = 0;
	
	/// Magnetic field at a certain point from a sub-domain
	virtual vec3 fieldAt(const vec3& v, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	/// Magnetic field with transform at a certain point from a sub-domain
	virtual vec2 fieldAtWithTransform(const vec3& v, const Matrix<2,3,mdouble>& M, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const { return fieldAt(v,vec2(0,0),vec2(1,1)); }
	/// Magnetic field at a point, with interaction matrix
	virtual vec2 fieldAtWithTransform(const vec3& v, const Matrix<2,3,mdouble>& M) const { return fieldAtWithTransform(v,M,vec2(0,0),vec2(1,1)); }
	
	/// Display field contributions over grid to given point
	void displayContribGrid(const vec3& v, unsigned int nx = 7, unsigned int ny = 7) const;
	
	unsigned int dflt_integrator_ndivs_x;	//< default number of sections to partition x integral in
	unsigned int dflt_integrator_ndivs_y;	//< default number of sections to partition y integral in
	
protected:
	mutable Integrator2D myIntegrator;	//< surface field integrator
};


#endif
