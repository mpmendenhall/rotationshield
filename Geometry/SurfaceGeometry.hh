#ifndef SURFACEGEOMETRY_HH
#define SURFACEGEOMETRY_HH 1

#include "Typedefs.hh"
#include "DVFunc.hh"
#include "Integrator.hh"

/// base class for defining surface mapping [0,1]^2 -> R^3
class SurfaceGeometry: public DVFunc<2,3,double> {
public:
	/// constructor
	SurfaceGeometry() { myIntegrator.setMethod(INTEG_GSL_QAG); }
	/// destructor
	virtual ~SurfaceGeometry() {}

	/// evaluate function
	virtual vec3 operator()(const vec2& x) const = 0;
	
	/// surface normal at point p; when not normalized, magnitude is dA
	virtual vec3 snorm(const vec2& l, bool normalized = false) const;
	
	/// pathlength for moves dl1, dl2 at surface position l
	virtual vec2 d_pathlength(vec2 l) const;
	
	/// differential area dA / dl1 dl2
	virtual double dA(const vec2& l) const;
	
	/// surface area corresponding to l range ll to ur
	virtual double area(const vec2& ll, const vec2& ur);
	
	/// calculate rotation to surface local coordinates v_local = M*v
	Matrix<3,3,double> rotToLocal(const vec2& x) const;
	
	/// calculate min and max distance to corners of region
	void proximity(vec3 x, vec2 ll, vec2 ur, double& mn, double& mx) const;
	
	/// whether surface is closed/periodic along particular axis
	virtual bool isClosed(unsigned int) const { return false; }
	
protected:

	Integrator myIntegrator;		//< integrator for internal calculations
	IntegratorND myIntegratorND;	//< multidimensional integrator for internal calculations
	
	// cache trig functions, since repeated calls will likely be for same values
	void cache_sincos(double theta, double& s, double& c) const;
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
	virtual vec3 deriv(const vec2& x, unsigned int i) const;
	
	/// differential area dA / dl1 dl2
	virtual double dA(const vec2& l) const;
	
	/// surface area corresponding to l range ll to ur
	virtual double area(const vec2& ll, const vec2& ur);
	
	// geometry-defining functions
	DVFunc1<2,double>* zr_profile;	//< z,r (l) profile
	
	/// whether surface is closed/periodic along particular axis
	virtual bool isClosed(unsigned int a) const { return a==0 ? (zr_profile && zr_profile->period) : a==1; }
	
protected:
	// cache profile calls, since likely to be for same value
	vec2 cache_profile(double l) const;
};


#endif
