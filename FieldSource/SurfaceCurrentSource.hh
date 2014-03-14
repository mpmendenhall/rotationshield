/// \file SurfaceCurrentSource.hh \brief Magnetic field from a surface current density

#ifndef SURFACECURRENTSOURCE_HH
/// Makes sure to only load this file once
#define SURFACECURRENTSOURCE_HH 1

#include "SurfaceSource.hh"
#include "SurfaceGeometry.hh"

class SurfaceCurrentSource: public SurfaceSource {
public:
	/// constructor
	SurfaceCurrentSource(SurfaceGeometry* SG = NULL): SurfaceSource(), mySurface(SG), vis_n1(200), vis_n2(200) {}
	
	/// current element dI from surface coordinate l
	vec3 dI_contrib(const vec2& l) const;
	
	/// field contribution f(x,y)dA; x,y in [0,1]^2
	virtual vec3 fieldAt_contrib_from(const vec3& v, const vec2& l) const;

	/// net current from specified region
	virtual vec3 netCurrent(vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
	/// Net current of the FieldSource
	virtual vec3 dipole() const { return netCurrent(vec2(0,0),vec2(1,1)); }
	
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
