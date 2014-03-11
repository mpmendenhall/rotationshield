/// \file SurfaceCurrentSource.hh \brief Magnetic field from a surface current density

#ifndef SURFACECURRENTSOURCE_HH
/// Makes sure to only load this file once
#define SURFACECURRENTSOURCE_HH 1

#include "SurfaceSource.hh"
#include "SurfaceGeometry.hh"

class SurfaceCurrentSource: public SurfaceSource {
public:
	/// constructor
	SurfaceCurrentSource(SurfaceGeometry* SG = NULL): SurfaceSource(), mySurface(SG), vis_n1(100), vis_n2(200) {}
	
	/// current element dipole dl from surface coordinate x,y
	vec3 dipoleContrib(mdouble x, mdouble y, vec3& xout) const;
	
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, mdouble x, mdouble y) const;

	SurfaceGeometry* mySurface;	//< surface over which current is distributed
	vec2 (*sj)(vec2, void*);	//< surface current density as a function of surface coordinate
	void* sjparams;				//< extra parameters for surface function
	
	/// Visualize the field source
	virtual void _visualize() const;
	unsigned int vis_n1;		//< visualization gridding, z
	unsigned int vis_n2;		//< visualization gridding, phi
	
	/// visualize local coordinate axes
	virtual void vis_coords(const vec2& l, double s = 0.02) const;
	
};


#endif
