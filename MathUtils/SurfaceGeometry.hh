#ifndef SURFACEGEOMETRY_HH
#define SURFACEGEOMETRY_HH 1

#include "Typedefs.hh"
#include "DVFunc.hh"

/// base class for defining surface mapping [0,1]^2 -> R^3
class SurfaceGeometry: public DVFunc<2,3,mdouble> {
public:
	/// constructor
	SurfaceGeometry() {}
	/// destructor
	virtual ~SurfaceGeometry() {}

	/// evaluate function
	virtual vec3 operator()(const vec2& x) const = 0;
	
	/// surface normal at point p; when not normalized, magnitude is dA
	virtual vec3 snorm(const vec2& p, bool normalized = false) const;
	
	/// calculate rotation to surface local coordinates v_local = M*v
	Matrix<3,3,mdouble> rotToLocal(const vec2& x) const;
};

/// cylindrically symmetric surface geometry
class CylSurfaceGeometry: public SurfaceGeometry {
public:
	/// constructor
	CylSurfaceGeometry(): fprofile(NULL) {}
	
	/// evaluate function
	virtual vec3 operator()(const vec2& x) const;
	
	/// partial derivative along axis i
	virtual vec3 deriv(const vec2& x, unsigned int i) const;
	
	// geometry-defining functions
	DVFunc1<2,mdouble>* fprofile;	//< z,r (l) profile
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
