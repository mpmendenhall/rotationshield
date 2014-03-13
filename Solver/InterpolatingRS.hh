/// \file InterpolatingRS.hh \brief Virtual base class for interpolated collections of interacting degrees of freedom

#ifndef INTERPOLATINGRS_HH
/// Makes sure to only load this file once
#define INTERPOLATINGRS_HH 1

#include "ReactiveSet.hh"
#include "InterpolationHelper.hh"
#include <vector>

/// ReactiveSet with utilities for interpolating between DF
class InterpolatingRS: public ReactiveSet {
public:
	/// constructor
	InterpolatingRS(unsigned int nph): ReactiveSet(nph) {}
		
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return InterplDF.n_pts(); }
	
protected:

	/// additional routines for setting a DF value
	virtual void _setDF(unsigned int DF, double v) { InterplDF[DF] = v; }
	/// additional routines for setting entire state vector
	virtual void _setDF(const mvec& v) { assert(v.size() == nDF()); InterplDF.setData(&v[0]); }
	
	InterpolationHelper InterplDF;	//< degrees of freedom in InterpolatingHelper grid
};

#endif
