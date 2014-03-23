#ifndef FIELDADAPTIVESURFACE_HH
#define FIELDADAPTIVESURFACE_HH 1

#include "DVFunc.hh"
#include "FieldEstimator2D.hh"
#include "BicubicGrid.hh"


/// Profile for optimizing cylinder to field
class FieldAdaptiveSurface: public DVFunc1<2,double> {
public:
	/// constructor
	FieldAdaptiveSurface(const DVFunc1<2,double>& f);
	
	/// evaluate function
	virtual vec2 operator()(double x) const { if(F.period) x=wrap(x); return F(l_remap(x)); }
	
	/// derivative
	virtual vec2 deriv(double x) const { if(F.period) x=wrap(x); return F.deriv(l_remap(x)) * l_remap.deriv(x); }
	
	/// optimize for field configuration
	void optimizeSpacing(const FieldEstimator2D& fes, double pfixed = 0.5, bool useDeriv = true);
	/// set constant spacing
	void setConstantSpacing();
	
	/// print out test points
	void symmetry_test() const;
	
protected:
	/// derivative of l distortion parameter
	double l_dist_deriv(double l) const;
	
	/// wrap a number into [0,1) for periodic function, starting half-way on outside
	double wrap(double x) const { double i; return x>=0 ? modf(x,&i) : 1+modf(x,&i); }
	
	CubicGrid l_remap;	//< distortion function, as interpolator
	
	const DVFunc1<2,double>& F;	//< reference function being distorted
};

#endif
