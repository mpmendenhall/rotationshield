#ifndef SURFACEGEOMETRY_HH
/// Make sure this header is only loaded once
#define SURFACEGEOMETRY_HH

#include "Typedefs.hh"
#include "DVFunc.hh"
#include "Integrator.hh"
#include <gsl/gsl_multimin.h>

/// base class for defining surface mapping [0,1]^2 -> R^3
class SurfaceGeometry: public DVFunc<2,3,double> {
public:
	/// constructor
	SurfaceGeometry();
	/// destructor
	virtual ~SurfaceGeometry();

	/// evaluate function
	virtual vec3 operator()(const vec2& x) const = 0;
	
	/// surface normal at point p; when not normalized, magnitude is dA
	virtual vec3 snorm(const vec2& l, bool normalized = false) const;
	
	/// pathlength for moves dl1, dl2 at surface position l
	virtual vec2 d_pathlength(vec2 l) const;
	
	/// differential area dA / dl1 dl2
	virtual double dA(const vec2& l) const;
	
	/// surface area corresponding to l range ll to ur
	virtual double area(const vec2& ll, const vec2& ur) const;
	
	/// calculate rotation to surface local coordinates v_local = M*v
	Matrix<3,3,double> rotToLocal(const vec2& x) const;
	
	/// calculate min and max distance to corners of region
	void proximity(vec3 x, vec2 ll, vec2 ur, double& mn, double& mx) const;
	
	/// find closest point on surface, and distance^2
	vec2 closestPoint(vec3 x, double& d2) const;
	
	
	/// whether surface is closed/periodic along particular axis
	virtual bool isClosed(unsigned int) const { return false; }
	
	unsigned int dflt_integrator_ndivs_x;	///< default number of sections to partition x integral in
	unsigned int dflt_integrator_ndivs_y;	///< default number of sections to partition y integral in
	vec2* polar_integral_center;			///< optional center point for switching to polar mode
	double polar_r0;						///< starting radius for polar integrals
	
	/// convenience mechanism for surface integrals
	mvec subdividedIntegral(mvec (*f)(vec2, void*), unsigned int fdim, void* fparams, vec2 ll, vec2 ur, unsigned int ndx=0, unsigned int ndy=0) const;

	mutable Integrator myIntegrator;		///< integrator for internal calculations
	IntegratorND myIntegratorND;			///< multidimensional integrator for internal calculations
	mutable Integrator2D myIntegrator2D;	///< surface field integrator, for polar-form CQUAD integrations

protected:
	
	// cache trig functions, since repeated calls will likely be for same values
	void cache_sincos(double theta, double& s, double& c) const;
	
	gsl_multimin_fdfminimizer* myMinimizer;	///< minimizer routine as needed
};

/// Convenience 2D vector surface
class Plane3D: public SurfaceGeometry {
public:
	Plane3D(vec3 lx = vec3(1,0,0), vec3 ly = vec3(0,1,0), vec3 o = vec3(-0.5,-0.5,0)): dx(lx), dy(ly), origin(o) {}
	virtual vec3 operator()(const vec2& l) const { return origin + dx*l[0] + dy*l[1]; }
	vec3 dx;
	vec3 dy;
	vec3 origin;
};



/// cylindrically symmetric surface geometry
class CylSurfaceGeometry: public SurfaceGeometry {
public:
	/// constructor
	CylSurfaceGeometry(DVFunc1<2,double>* f = NULL): zr_profile(f) {}
	
	/// evaluate function
	virtual vec3 operator()(const vec2& x) const;
	
	/// partial derivative along axis i
	virtual vec3 deriv(const vec2& x, unsigned int i, bool normalized = false) const;
	
	/// differential area dA / dl1 dl2
	virtual double dA(const vec2& l) const;
	
	/// surface normal at point p; when not normalized, magnitude is dA
	virtual vec3 snorm(const vec2& l, bool normalized = false) const;
	
	/// surface area corresponding to l range ll to ur
	virtual double area(const vec2& ll, const vec2& ur) const;
	
	// geometry-defining functions
	DVFunc1<2,double>* zr_profile;	///< z,r (l) profile
	
	/// whether surface is closed/periodic along particular axis
	virtual bool isClosed(unsigned int a) const { return a==0 ? (zr_profile && zr_profile->period) : a==1; }
	
protected:
	// cache profile calls, since likely to be for same value
	vec2 cache_profile(double l) const;
};


/// Geometry with translation and matrix transformation
class TransformedGeometry: public SurfaceGeometry {
public:
	/// constructor
	TransformedGeometry(SurfaceGeometry* G): transf_M(Matrix<3,3,double>::identity()), baseGeom(G) {}
	
	/// evaluate function
	virtual vec3 operator()(const vec2& x) const { return transf_M * (*baseGeom)(x) + transl_v; }
	
	/// partial derivative along axis i
	virtual vec3 deriv(const vec2& x, unsigned int i, bool normalized = false) const { return (transf_M * baseGeom->deriv(x,i,normalized)).normalized(); }
	
	/// differential area dA / dl1 dl2
	virtual double dA(const vec2& l) const { return baseGeom->dA(l); } // TODO do this right
	
	/// whether surface is closed/periodic along particular axis
	virtual bool isClosed(unsigned int a) const { return baseGeom->isClosed(a); }
	
	vec3 transl_v;					///< translation vector
	Matrix<3,3,double> transf_M;	///< transformation matrix
	
protected:
	SurfaceGeometry* baseGeom;
};

#endif
