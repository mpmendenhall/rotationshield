#ifndef UNIFORMFIELD_HH
/// Makes sure to only load this file once
#define UNIFORMFIELD_HH 1

#include "FieldSource.hh"
#include <stdio.h>

/// A uniform magnetic field
class UniformField: public FieldSource {
public:
	/// Constructor
	UniformField(const vec3& f = vec3()): FieldSource(), B(f) {}
	/// Destructor
	virtual ~UniformField() {}
	
	vec3 B; //< the uniform field vector
	
	/// Magnetic field at a specified point
	virtual vec3 fieldAt(const vec3& v) const { return B; }
	/// Magnetic field averaged over specified line
	virtual vec3 fieldOverLine(Line l) const { return B; }
	/// Magnetic field averaged over specified plane
	virtual vec3 fieldOverPlane(Plane pl) const { return B; }

	/// Print info to stdout
	virtual void display() const { printf("[Uniform Field]\n"); }
};

#endif
