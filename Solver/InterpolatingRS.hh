/// \file InterpolatingRS.hh \brief Virtual base class for interpolated collections of interacting degrees of freedom

#ifndef INTERPOLATINGRS_HH
/// Makes sure to only load this file once
#define INTERPOLATINGRS_HH 1

#include "ReactiveSet.hh"
#include "InterpolationHelper.hh"
#include "BicubicGrid.hh"
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
	virtual void _setDFv(const mvec& v) { assert(v.size() == nDF()); InterplDF.setData(&v[0]); }
	
	InterpolationHelper InterplDF;	//< degrees of freedom in InterpolatingHelper grid
};

/// Simplified case for 2D interpolation
class InterpolatingRS2D: public ReactiveSet {
public:
	/// constructor
	InterpolatingRS2D(unsigned int nph): ReactiveSet(nph), nZ(0), nDFi(0) { }
	/// destructor
	virtual ~InterpolatingRS2D() { clear_data(); }
	
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return nZ*nPhi*nDFi; }
	
protected:
	/// determine "address" for DF in interpolating grids
	void DF_address(unsigned int DF, unsigned int& p, unsigned int& z, unsigned int& d) const;
	
	/// additional routines for setting a DF value
	virtual void _setDF(unsigned int DF, double v);
	
	//virtual void _setDFv(const mvec& v) { ReactiveSet::_setDFv(v); G[0]->printData(); }
	
	/// delete previous grids
	virtual void clear_data();
	/// set up data grid
	virtual void make_grids(unsigned int nz, unsigned int ndf);
	
	std::vector<BicubicGrid*> G;	//< degrees of freedom stored in interpolator grid
	unsigned int nZ;				//< grid size in z direction
	unsigned int nDFi;				//< number of DF per element
};


#endif
