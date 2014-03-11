/// \file SurfaceCurrentRS.hh \brief ReactiveSet for a surface current responding to magnetic field conditions

#ifndef SURFACECURRENTRS_HH
/// Makes sure to only load this file once
#define SURFACECURRENTRS_HH 1

#include "SurfaceCurrentSource.hh"
#include "InterpolatingRS.hh"

/// defines surface boundary condition by producing surface current in response to applied magnetic field
class SurfaceI_Response {
public:
	/// constructor
	SurfaceI_Response(mdouble mu = 0): murel(mu) { setRmat(); }
	
protected:
	mdouble murel;	//< relative permeability
	mmat rmat;		//< response matrix to applied field
	
	/// generate correct response matrix to applied fields
	void setRmat() {
		//Matrix<2,3,mdouble> M = Matrix<2,3,mdouble>();
		//M(0,1) = M(1,0) = 2.0*(1.0-murel)/(1.0+murel);
		//rmat = mmat(M*p.projectionMatrix());
	}
};


/// Continuous surface current responding to magnetic field
/// things needing to be set to implement:
///		SurfaceGeometry* mySurface
class SurfaceCurrentRS: public SurfaceCurrentSource, public InterpolatingRS {
public:
	/// constructor
	SurfaceCurrentRS(unsigned int nph, unsigned int nz);

	/// set interaction protocol to use; respond whether accepted
	virtual bool set_protocol(void* ip);
	/// respond to interaction protocol
	virtual void queryInteraction();
	
	/// set surface response at all points
	void setSurfaceResponse(SurfaceI_Response r);
	
	/// get interpolated surface response
	vec2 eval(const vec2& p) const;
	
	const unsigned int nZ;	//< number of z divisions
	
	/// visualization routine
	virtual void _visualize() const;
	
protected:
	std::vector<SurfaceI_Response> sdefs;
};

#endif