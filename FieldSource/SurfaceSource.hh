/// \file SurfaceSource.hh \brief Base class for extended 2D surface magnetic field sources

#ifndef SURFACESOURCE_HH
/// Makes sure to only load this file once
#define SURFACESOURCE_HH 1

#include "FieldSource.hh"
#include "Integrator.hh"
#include "SurfaceGeometry.hh"

class SurfaceSource: public FieldSource {
public:
	/// constructor
	SurfaceSource(SurfaceGeometry* SG = NULL, const std::string& nm = "SurfaceSource"):
		FieldSource(nm), mySurface(SG), dflt_integrator_ndivs_x(1), dflt_integrator_ndivs_y(1), polar_integral_center(NULL), polar_r0(0) {}

	/// destructor
	virtual ~SurfaceSource() {}
	
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, const vec2& l) const = 0;
	
	/// Magnetic field at a certain point from a sub-domain
	virtual vec3 fieldAt(const vec3& v, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	/// Magnetic field with transform at a certain point from a sub-domain
	virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,mdouble>& M, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	/// Magnetic field with transform at a certain point from a sub-domain
	virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,mdouble>& M, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const { return fieldAt(v,vec2(0,0),vec2(1,1)); }
	/// Magnetic field at a point, with interaction matrix
	virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,mdouble>& M) const { return fieldAtWithTransform2(v,M,vec2(0,0),vec2(1,1)); }
	/// Magnetic field at a point, with interaction matrix
	virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,mdouble>& M) const { return fieldAtWithTransform3(v,M,vec2(0,0),vec2(1,1)); }
		
	/// Display field contributions over grid to given point
	void displayContribGrid(const vec3& v, unsigned int nx = 7, unsigned int ny = 7) const;
	
	SurfaceGeometry* mySurface;	//< surface on which source is defined


	unsigned int dflt_integrator_ndivs_x;	//< default number of sections to partition x integral in
	unsigned int dflt_integrator_ndivs_y;	//< default number of sections to partition y integral in
	
protected:

	mutable Integrator2D myIntegrator;	//< surface field integrator
	vec2* polar_integral_center;		//< optional center point for switching to polar mode
	double polar_r0;					//< starting radius for polar integrals
	
	/// convenience mechanism for integrations split over may divisions, limited to surface [0,1]->[0,1] // TODO the problem??
	mvec subdividedIntegral(mvec (*f)(vec2, void*), void* fparams, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const;
};


#endif
