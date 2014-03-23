#ifndef TYPEDEFS_HH
#define TYPEDEFS_HH 1

#include "Vec.hh"
#include "Matrix.hh"
#include "VarVec.hh"
#include "VarMat.hh"

typedef Vec<4,double> vec4;
typedef Vec<3,double> vec3;
typedef Vec<2,double> vec2;
typedef Matrix<3,3,double> mat3;

typedef VarMat<double> mmat;
typedef VarVec<double> mvec;

template<typename T>
inline T sign(T t) { if(t<0) return -1; return (t==0)?0:1; }

#endif
