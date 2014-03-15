/// \file FieldSource.hh \brief Base class for magnetic field sources

#ifndef FIELDSOURCE_HH
/// Makes sure to only load this file once
#define FIELDSOURCE_HH 1

#include "Geometry.hh"
#include "Typedefs.hh"
#include "RefCounter.hh"
#include "Integrator.hh"
#include "Matrix.hh"

/// Base (virtual) class for magnetic field sources due to current distributions
class FieldSource: public RefCounter, public Visualizable {
public:
	/// Constructor
	FieldSource(): RefCounter() {}
	/// Destructor
	virtual ~FieldSource() {}
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const = 0;
	/// Magnetic field at a point, with interaction matrix (optionally subclass for improved numerics)
	virtual vec2 fieldAtWithTransform(const vec3& v, const Matrix<2,3,mdouble>& M) const { return M*fieldAt(v); }
	/// Magnetic field averaged over specified line
	virtual vec3 fieldOverLine(Line l) const;
	/// Magnetic field averaged over specified plane
	virtual vec3 fieldOverPlane(Plane pl) const;
	/// Calculate effective net current of field source
	virtual vec3 net_current() const { return vec3(); }
	/// Print info to stdout
	virtual void display() const { printf("[FieldSource]\n"); }
		
private:
	static Integrator lineIntegrator; //< Integrator for averaging over lines
	static Integrator planeIntegrator; //< Integrator for averaging over planes
};

#endif
