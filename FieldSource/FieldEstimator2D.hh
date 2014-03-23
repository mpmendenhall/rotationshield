/// \file FieldEstimator2D.hh \brief 2D magnetic fields

#ifndef FIELDESTIMATOR2D_HH
/// Makes sure to only load this file once
#define FIELDESTIMATOR2D_HH 1

#include "FieldSource.hh"
#include <vector>

/// Rough estimate of field magnitude (in 2D plane) due to wires perpendicular to shield, used to optimize shield gridding for cos theta coil endcaps
class FieldEstimator2D {
public:
	/// Constructor
	FieldEstimator2D() {}
	/// Destructor
	virtual ~FieldEstimator2D() {}
	/// Estimated field at a point
	virtual vec2 estimateAt(const vec2& v) const;
	/// estimate rate of change in given direction
	virtual vec2 derivAt(const vec2& v, const vec2& dx) const;
	/// Add a "line source" perpendicular to plane of interest
	void addsource(const vec2& v, double j) { sources.push_back(v); currents.push_back(j); }
protected:
	std::vector<vec2> sources;
	std::vector<double> currents;
};

/// Field estimator based off 3D field source
class FieldEstimator2Dfrom3D: public FieldEstimator2D {
public:
	/// Constructor
	FieldEstimator2Dfrom3D(FieldSource* FS): myFS(FS) { assert(myFS); myFS->retain(); }
	/// Destructor
	virtual ~FieldEstimator2Dfrom3D() { myFS->release(); }
	/// Estimated field at a point
	virtual vec2 estimateAt(const vec2& v) const;
protected:
	FieldSource* myFS;
};

#endif
