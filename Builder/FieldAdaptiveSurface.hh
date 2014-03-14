#ifndef FIELDADAPTIVESURFACE_HH
#define FIELDADAPTIVESURFACE_HH 1

#include "DVFunc.hh"
#include "FieldEstimator2D.hh"
#include "BicubicGrid.hh"


/// Profile for optimizing cylinder to field
class FieldAdaptiveSurface: public DVFunc1<2,mdouble> {
public:
	/// constructor
	FieldAdaptiveSurface(const DVFunc1<2,mdouble>& f);
	
	/// evaluate function
	virtual vec2 operator()(mdouble x) const { return F(l_remap(x)); }
	
	/// derivative
	virtual vec2 deriv(mdouble x) const { return F.deriv(l_remap(x)) * l_remap.deriv(x); }
	
	/// optimize for field configuration
	void optimizeSpacing(const FieldEstimator2D& fes, double pfixed = 0.5, bool useDeriv = true);
	/// set constant spacing
	void setConstantSpacing();
	
	/// print out test points
	void symmetry_test() const;
	
protected:
	/// derivative of l distortion parameter
	mdouble l_dist_deriv(mdouble l) const;
	
	CubicGrid l_remap;	//< distortion function, as interpolator
	
	const DVFunc1<2,mdouble>& F;	//< reference function being distorted
};

#endif
