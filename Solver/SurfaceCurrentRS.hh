/// \file SurfaceCurrentRS.hh \brief ReactiveSet for a surface current responding to magnetic field conditions

#ifndef SURFACECURRENTRS_HH
/// Makes sure to only load this file once
#define SURFACECURRENTRS_HH 1

#include "SurfaceCurrentSource.hh"
#include "InterpolatingRS.hh"
#include "MagRS.hh"

/// defines surface boundary condition by producing surface current in response to applied magnetic field
class SurfaceI_Response {
public:
	/// constructor
	SurfaceI_Response(mdouble mu = 0) { setMu(mu); }
	
	/// generate correct response matrix to applied fields for given relative permeability
	void setMu(double mu) {
		murel = mu;
		rmat2 = Matrix<2,3,mdouble>();
		rmat2(0,1) = 2*(murel-1)/(murel+1);
		rmat2(1,0) = -rmat2(0,1);
	}
	
	mdouble murel;					//< relative permeability
	Matrix<2,3,mdouble> rmat2;		//< response matrix to applied field
};


/// Continuous surface current responding to magnetic field
class SurfaceCurrentRS: public MagF_Responder, public SurfaceCurrentSource, public InterpolatingRS2D {
public:
	/// constructor; set xdf=1 to enable dipole density response
	SurfaceCurrentRS(unsigned int nph, unsigned int nz, const std::string& nm = "SurfaceCurrentRS");

	// ReactiveSet subclassed functions
	//=====================================
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0);
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip);
	//=====================================
	
	/// calculate response to incident field
	virtual void calculateIncident(const FieldSource& f);
	
	/// set surface response at all points
	void setSurfaceResponse(SurfaceI_Response r);
	
	/// visualization routine
	virtual void _visualize() const;
	/// visualize current vectors
	virtual void vis_i_vectors(double s = 0.5, double mx = 0.2) const;
	
	/// get center surface coordinates for i^th element
	vec2 surf_coords(unsigned int i) const { return vec2( ((i/nPhi)+0.5)/nZ, ((i%nPhi)+0.5)/nPhi); }
	/// get surface coordinate range for element (not including extended interpolation effects)
	void element_surface_range(unsigned int i, vec2& ll, vec2& ur) const;
	
	/// get interpolated surface current response
	vec2 eval_J(const vec2& p) const { return vec2( (*G[0])(p[0],p[1]), (*G[1])(p[0],p[1]) ); }	
	
protected:
	
	/// one surface element's reaction
	mvec subelReaction(ReactiveSet* R);
	
	unsigned int ixn_el;	//< interacting element currently being probed
	
	std::vector<SurfaceI_Response> sdefs;	//< surface response definitions at each site
};

#endif