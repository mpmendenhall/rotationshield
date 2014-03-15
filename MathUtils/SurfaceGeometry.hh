#ifndef SURFACEGEOMETRY_HH
#define SURFACEGEOMETRY_HH 1

#include "Typedefs.hh"
#include "DVFunc.hh"
#include "Integrator.hh"

/// base class for defining surface mapping [0,1]^2 -> R^3
class SurfaceGeometry: public DVFunc<2,3,mdouble> {
public:
	/// constructor
	SurfaceGeometry() { myIntegrator.setMethod(INTEG_GSL_QAG); }
	/// destructor
	virtual ~SurfaceGeometry() {}

	/// evaluate function
	virtual vec3 operator()(const vec2& x) const = 0;
	
	/// surface normal at point p; when not normalized, magnitude is dA
	virtual vec3 snorm(const vec2& l, bool normalized = false) const;
	
	/// differential area dA / dl1 dl2
	virtual mdouble dA(const vec2& l) const;
	
	/// surface area corresponding to l range ll to ur
	virtual mdouble area(const vec2& ll, const vec2& ur);
	
	/// calculate rotation to surface local coordinates v_local = M*v
	Matrix<3,3,mdouble> rotToLocal(const vec2& x) const;
	
protected:

	Integrator2D myIntegrator;	//< integrator for internal calculations
	
	// cache trig functions, since repeated calls will likely be for same values
	void cache_sincos(double theta, double& s, double& c) const;
};

/// cylindrically symmetric surface geometry
class CylSurfaceGeometry: public SurfaceGeometry {
public:
	/// constructor
	CylSurfaceGeometry(DVFunc1<2,mdouble>* f = NULL): zr_profile(f) {}
	
	/// evaluate function
	virtual vec3 operator()(const vec2& x) const;
	
	/// partial derivative along axis i
	virtual vec3 deriv(const vec2& x, unsigned int i) const;
	
	/// differential area dA / dl1 dl2
	virtual mdouble dA(const vec2& l) const;
	
	/// surface area corresponding to l range ll to ur
	virtual mdouble area(const vec2& ll, const vec2& ur);
	
	// geometry-defining functions
	DVFunc1<2,mdouble>* zr_profile;	//< z,r (l) profile
	
protected:
	// cache profile calls, since likely to be for same value
	vec2 cache_profile(mdouble l) const;
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

/// Convenience linear sweep DFunc
class Line2D: public DVFunc1<2,mdouble> {
public:
	Line2D(vec2 a, vec2 b): x0(a), x1(b) {}
	vec2 x0,x1;
	virtual vec2 operator()(mdouble x) const { return x0*(1-x) + x1*x; }
	virtual vec2 deriv(mdouble) const { return x1-x0; }
};

#endif
