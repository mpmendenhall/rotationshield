#ifndef FIELDADAPTIVESURFACE_HH
#define FIELDADAPTIVESURFACE_HH 1

#include "DVFunc.hh"
#include "FieldEstimator2D.hh"
#include "InterpolationHelper.hh"


/// Profile for optimizing cylinder to field
class FieldAdaptiveSurface: public DVFunc1<2,mdouble> {
public:
	/// constructor
	FieldAdaptiveSurface(const DVFunc1<2,mdouble>& f);
	
	/// evaluate function
	virtual vec2 operator()(mdouble x) const;
	
	/// derivative
	virtual vec2 deriv(mdouble x) const;
	
	/// optimize for field configuration
	void optimizeSpacing(const FieldEstimator2D& fes, double pfixed = 0.5);
	/// set constant spacing
	void setConstantSpacing();
	
protected:
	/// derivative of l distortion parameter
	mdouble l_dist_deriv(mdouble l) const;
	
	InterpolationHelper Ilz;	//< distortion function, as interpolator
	
	const DVFunc1<2,mdouble>& F;	//< reference function being distorted
};

#endif
