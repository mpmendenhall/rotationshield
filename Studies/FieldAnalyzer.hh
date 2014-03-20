/// \file "FieldAnalyzer.hh" \brief Routines for analyzing the simulated fields
#ifndef ANALYSIS_HH
/// Make sure this header is only loaded once
#define ANALYSIS_HH 1

#include "FieldSource.hh"
#include "VarVec.hh"
#include "linmin.hh"
#include "Geometry.hh"
#include <iostream>
#include <iomanip>

/// A class for analyzing field uniformity/shape over rectangular grids.
///
/// Uses Simpson's Rule to numerically integrate the fields
///	over the specified volume --- good enough for the smoothly
///	varying fields typically encountered.
class FieldAnalyzer {
public:
	/// Constructor \param fs the FieldSource to analyze fields from
	FieldAnalyzer(FieldSource* fs): FS(fs) { FS->retain(); }
	/// destructor
	~FieldAnalyzer() { FS->release(); }
	/// Survey the fields over a rectangular grid
	/// \param ll lower-left corner of the survey area
	/// \param ur epur-right corner of the survey area
	/// \param nX number of grid points along the x direction
	/// \param nY number of grid points along the y direction
	/// \param nZ number of grid points along the z direction
	/// \param datf file to store surveyed fields to
	/// \param statsout ostream to which summary fields data is written
	/// \param key string to prepend to summary data
	void survey(vec3 ll, vec3 ur, unsigned int nX, unsigned int nY, unsigned int nZ, std::ostream& statsout = std::cout, std::ostream& datsout = std::cout) const;
	
	void visualizeSurvey(vec3 ll, vec3 ur, unsigned int nX, unsigned int nY, unsigned int nZ) const;
	
	/// Get the Simpson's Rule numerical integrating weight for a point
	static mdouble simpsonCoeff(unsigned int n, unsigned int ntot);

protected:
	const FieldSource* FS; //< the FieldSource producing the fields being analyzed
};

#endif
