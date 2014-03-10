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
	
	/// surface normal at point p; magnitude is dA / dl1 dl2
	virtual vec3 snorm(const vec2& p) const;
	

};

/// cylindrically symmetric surface geometry
class CylSurfaceGeometry: public SurfaceGeometry {
public:
	/// constructor
	CylSurfaceGeometry(): fz(NULL), fr(NULL) {}
	
	/// evaluate function
	virtual vec3 operator()(const vec2& x) const;
	
	/// partial derivative along axis i
	virtual vec3 deriv(const vec2& x, unsigned int i) const;
	
	// geometry-defining functions
	DFunc<mdouble>* fz;	//< z position as a function of l1 in [0,1]
	DFunc<mdouble>* fr;	//< radius as a function of z position
};


#endif
