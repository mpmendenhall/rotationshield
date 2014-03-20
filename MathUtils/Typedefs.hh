#ifndef TYPEDEFS_HH
#define TYPEDEFS_HH 1

#include "Vec.hh"
#include "Matrix.hh"
#include "VarVec.hh"
#include "VarMat.hh"

typedef double mdouble; //< Sets precision at which most calculations are carried out

typedef Vec<4,mdouble> vec4;
typedef Vec<3,mdouble> vec3;
typedef Vec<2,mdouble> vec2;
typedef Matrix<3,3,mdouble> mat3;

typedef VarMat<mdouble> mmat;
typedef VarVec<mdouble> mvec;

template<typename T>
inline T sign(T t) { if(t<0) return -1; return (t==0)?0:1; }

#endif
