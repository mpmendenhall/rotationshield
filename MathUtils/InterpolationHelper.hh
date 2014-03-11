#ifndef INTERPOLATIONHELPER_HH
#define INTERPOLATIONHELPER_HH 1

#include "Interpolator.hh"

/// Multi-dimensional interpolating grid; also allows "flat" access to data points
class InterpolationHelper: protected DataSequence {
public:
	/// constructor
	InterpolationHelper();
	/// destructor
	virtual ~InterpolationHelper();
	
	/// set up n1-n0 - dimensional data grid with (*n0), ... points in each dimension
	void setupDataGrid(const unsigned int* n0, const unsigned int* n1);
	/// convenience form using std::vector
	void setupDataGrid(const std::vector<unsigned int>& n) { setupDataGrid(&*n.begin(), &*n.end()); }
	
	
	/// direct access to point by multi-dimensional index
	double& operator[](const unsigned int* n);
	/// direct access to point by flat index
	double& operator[](unsigned int n);
	/// zero out all data points
	void zeroData();
	/// set "flat" data array
	void setData(const double* x0);
	
	
	/// evaluation at given x
	double eval(const double* x) const;
	
	/// directly address sub-InterpolationHelper
	const InterpolationHelper& getSubHelper(unsigned int nDeep, const unsigned int* n) const;
	/// get a list of sub-interpolators at given depth
	std::vector<InterpolationHelper*> getSubHelpers(unsigned int nDeep);
	/// get interpolator
	Interpolator* getInterpolator() { return myInterpolator; }
	
	/// set interpolation method for all InterpolationHelpers at given depth
	void setInterpolatorMethod(Interpolator* (*makeInterp)(DataSequence*), unsigned int nDeep = 0);
	/// set boundary condition at indicated depth
	void setBoundaryCondition(BoundaryCondition b, unsigned int nDeep = 0);
	
	/// get number of underlying datapoints
	unsigned int n_pts() const { return cum_sum_dpts.back(); }
	/// recalculate point count structure
	void recalc_structure(bool recurse=true);

protected:

	/// value at given grid location
	virtual double valueAt(int i, void* xopts) const;
	
	/// delete sub-interpolators
	void clearSubinterps();
	
	/// check whether this is the bottom layer of interpolators
	bool isBottom() const { return !subInterpolators.size(); }
	
	/// calculate sub-index of "flat" index; return sub-interpolator #; replace n with its flat index
	unsigned int from_flat_index(unsigned int& n) const;
	
	Interpolator* myInterpolator;
	std::vector<InterpolationHelper*> subInterpolators;	//< sub-interpolators if not lowest dimension
	std::vector<double> myData;							//< internal list of data points if lowest dimension
	std::vector<unsigned int> cum_sum_dpts;				//< cumulative sum of datapoints in underlying structures
};


#endif
