/// \file FieldSource.hh \brief Base class for magnetic field sources

#ifndef FIELDSOURCE_HH
/// Makes sure to only load this file once
#define FIELDSOURCE_HH 1

#include "geometry.hh"
#include "MiscUtils.hh"
#include "RefCounter.hh"
#include "integrator.hh"
					
/// Base (virtual) class for magnetic field sources due to current distributions
class FieldSource: public RefCounter {
public:
	/// Constructor
	FieldSource(): RefCounter() {}
	/// Destructor
	virtual ~FieldSource() {}
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const = 0;
	/// Magnetic field averaged over specified line
	virtual vec3 fieldOverLine(Line l) const;
	/// Magnetic field averaged over specified plane
	virtual vec3 fieldOverPlane(Plane pl) const;
	/// Net current of the FieldSource
	virtual vec3 dipole() const { return vec3(); }
	/// Print info to stdout
	virtual void display() const { printf("[FieldSource]\n"); }
	
	/// Visualize the field source
	virtual void visualize(bool top = true, mdouble scale = 1.0) const { printf("Unimplemented!!!\n"); }
	
private:
	static Integrator lineIntegrator; //< Integrator for averaging over lines
	static Integrator planeIntegrator; //< Integrator for averaging over planes
};

#endif
