/* 
 * Typedefs.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

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
