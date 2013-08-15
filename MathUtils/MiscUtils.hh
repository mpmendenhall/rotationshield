/// \file "MiscUtils.hh" \brief Assorted useful functions
#ifndef MISCUTILS_HH
/// Make sure this header is only loaded once
#define MISCUTILS_HH 1

#include <cassert>
#include <string>
#include "Vec.hh"
#include "Matrix.hh"

/// numerical approximation for \f$ \pi \f$
#define PI 3.14159265358979323846264338327950288419716939937510

typedef double mdouble; //< Sets precision at which most calculations are carried out

template<typename T>
inline T max(T a, T b) { return (a>=b)?a:b; }
template<typename T>
inline T min(T a, T b) { return (a<=b)?a:b; }

template<typename T>
inline T sign(T t) { if(t<0) return -1; return (t==0)?0:1; }
	
inline float max(float a, float b) { return (b<a)?a:b; }
inline float min(float a, float b) { return (a<b)?a:b; }

/// normalize an angle to [t0,t0+2*PI)
float normalizeAngle(float a, float theta0 = -PI);

/// Return a random number, uniformly distributed over interval [a,b]
/**	\param a lower bound of interval
 \param b upper bound of interval
 \return a random number in the interval [a,b] */
mdouble randunif(mdouble a, mdouble b);

typedef Vec<4,mdouble> vec4;
typedef Vec<3,mdouble> vec3;
typedef Vec<2,mdouble> vec2;
typedef Matrix<3,3,mdouble> mat3;

#endif
