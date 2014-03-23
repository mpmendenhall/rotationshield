#include "InfiniteLineSource.hh"
#include <math.h>

vec3 InfiniteLineSource::fieldAt(const vec3& v) const {
	vec3 b = cross(l.s-v,l.dv);
	double m2 = b.mag2();
	if(!m2) return b;
	return b*(j/(4.0*M_PI*m2));
}
