/* 
 * Vec.hh, part of the RotationShield program
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

/// \file Vec.hh \brief Templatized fixed-length array class with mathematical operations
#ifndef VEC_HH
/// Make sure this header is only loaded once
#define VEC_HH

#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <iostream>
#include <vector>

/// Fixed-length vector arithmetic class
template<unsigned int N, typename T>
class Vec {
public:
	
	/// default constructor for zero vector
	Vec() { for(unsigned int i=0; i<N; i++) x[i] = T(); }
	/// constructor for 1-vectors
	Vec(T x0) { assert(N==1); x[0] = x0; }
	/// constructor for 2-vectors
	Vec(T x0, T x1) { assert(N==2); x[0] = x0; x[1] = x1; }
	/// constructor for 3-vectors
	Vec(T x0, T x1, T x2) { assert(N==3); x[0] = x0; x[1] = x1; x[2] = x2; }
	/// constructor for 4-vectors
	Vec(T x0, T x1, T x2, T x3) { assert(N==4); x[0] = x0; x[1] = x1; x[2] = x2; x[3] = x3; }
	/// constructor for 7-vectors (for, e.g., SI unit system)
	Vec(T x0, T x1, T x2, T x3, T x4, T x5, T x6 ) {
		assert(N==7);
		x[0] = x0; x[1] = x1; x[2] = x2;
		x[3] = x3; x[4] = x4; x[5] = x5;
		x[6] = x6;
	}
	
	static Vec<N,T> basis(unsigned int n) { Vec<N,T> v = Vec<N,T>(); v[n] = 1; return v; }
	
	/// print to stdout
	void display(const char* suffix = "\n") const;
	
	/// dot product with another vector
	T dot(const Vec<N,T>& v) const { T s = x[0]*v[0]; for(unsigned int i=1; i<N; i++) s+=x[i]*v[i]; return s; }
	/// square magnitude \f$ v \cdot v \f$
	T mag2() const { return dot(*this); }
	/// magnitude \f$ \sqrt{v\cdot v} \f$
	T mag() const { return sqrt(mag2()); }
	/// sum of vector elements
	T sum() const { T s = x[0]; for(unsigned int i=1; i<N; i++) s += x[i]; return s; }
	/// product of vector elements
	T prod() const { T s = x[0]; for(unsigned int i=1; i<N; i++) s *= x[i]; return s; }
	
	/// this vector, normalized to magnitude 1
	Vec<N,T> normalized() const { return (*this)/mag(); }
	/// project out component parallel to another vector
	Vec<N,T> paraProj(const Vec<N,T>& v) const { return v*(dot(v)/v.mag2()); }
	/// project out component orthogonal to another vector
	Vec<N,T> orthoProj(const Vec<N,T>& v) const { return (*this)-paraProj(v); }
	/// angle with another vector
	T angle(const Vec<N,T> v) const { return acos(dot(v)/sqrt(mag2()*v.mag2())); }
	
	/// mutable element access operator
	T& operator[](unsigned int i) { assert(i<N); return x[i]; }
	/// immutable element access operator
	const T& operator[](unsigned int i) const { assert(i<N); return x[i]; }
	
	/// unary minus operator
	const Vec<N,T> operator-() const;
	
	/// inplace addition
	Vec<N,T>& operator+=(const Vec<N,T>& rhs);	
	/// inplace addition of a constant
	Vec<N,T>& operator+=(const T& c);	
	/// inplace subtraction
	Vec<N,T>& operator-=(const Vec<N,T>& rhs);	
	/// inplace subtraction of a constant
	Vec<N,T>& operator-=(const T& c);
	
	/// inplace multiplication
	Vec<N,T>& operator*=(T c);
	/// inplace elementwise multiplication
	Vec<N,T>& operator*=(const Vec<N,T>& other);
	/// inplace division
	Vec<N,T>& operator/=(T c);
	/// inplace elementwise division
	Vec<N,T>& operator/=(const Vec<N,T>& other);
	
	/// addition operator
	const Vec<N,T> operator+(const Vec<N,T>& other) const;	
	/// subtraction operator
	const Vec<N,T> operator-(const Vec<N,T>& other) const;	
	
	/// multiplication operator
	const Vec<N,T> operator*(T c) const;
	/// elementwise multiplication operator
	const Vec<N,T> operator*(const Vec<N,T>& other) const;
	/// division operator
	const Vec<N,T> operator/(T c) const;
	/// elementwise division operator
	const Vec<N,T> operator/(const Vec<N,T>& other) const;
	
	/// equality operator
	bool operator==(const Vec<N,T>& rhs) const;
	/// inequality operator
	bool operator!=(const Vec<N,T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator<(const Vec<N,T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator<=(const Vec<N,T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator>(const Vec<N,T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator>=(const Vec<N,T>& rhs) const;
	
	/// write in binray form to a file
	void writeBinary(std::ostream& o) const { o.write((char*)x,N*sizeof(T)); }
	/// read a Vec from a file
	static Vec<N,T> readBinary(std::istream& s) { Vec<N,T> v; v.loadBinaryData(s); return v; }
	/// read in binary form from a file
	Vec<N,T>& loadBinaryData(std::istream& s) { s.read((char*)x,N*sizeof(T)); return *this; }

	
protected:
	T x[N]; ///< vector elements
};

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator-() const {
	Vec<N,T> v = *this;
	for(unsigned int i=0; i<N; i++)
		v[i] *= -1.0;
	return v;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator+=(const Vec<N,T>& rhs) {
	for(unsigned int i=0; i<N; i++)
		x[i] += rhs[i];
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator-=(const Vec<N,T>& rhs) {
	for(unsigned int i=0; i<N; i++)
		x[i] -= rhs[i];
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator+=(const T& c) {
	for(unsigned int i=0; i<N; i++)
		x[i] += c;
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator-=(const T& c) {
	for(unsigned int i=0; i<N; i++)
		x[i] -= c;
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator*=(T c) {
	for(unsigned int i=0; i<N; i++)
		x[i] *= c;
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator*=(const Vec<N,T>& other) {
	for(unsigned int i=0; i<N; i++)
		x[i] *= other[i];
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator/=(T c) {
	for(unsigned int i=0; i<N; i++)
		x[i] /= c;
	return *this;
}

template<unsigned int N, typename T>
Vec<N,T>& Vec<N,T>::operator/=(const Vec<N,T>& other) {
	for(unsigned int i=0; i<N; i++)
		x[i] /= other[i];
	return *this;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator+(const Vec<N,T>& other) const {
	Vec<N,T> result = *this;
	result += other;
	return result;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator-(const Vec<N,T>& other) const {
	Vec<N,T> result = *this;
	result -= other;
	return result;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator*(T c) const {
	Vec<N,T> result = *this;
	result *= c;
	return result;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator*(const Vec<N,T>& other) const {
	Vec<N,T> result = *this;
	result *= other;
	return result;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator/(T c) const {
	Vec<N,T> result = *this;
	result /= c;
	return result;
}

template<unsigned int N, typename T>
const Vec<N,T> Vec<N,T>::operator/(const Vec<N,T>& other) const {
	Vec<N,T> result = *this;
	result /= other;
	return result;
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator==(const Vec<N,T>& rhs) const {
	for(unsigned int i=0; i<N; i++)
		if(x[i] != rhs.x[i])
			return false;
	return true;
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator!=(const Vec<N,T>& rhs) const {
	return !(*this == rhs);
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator<(const Vec<N,T>& rhs) const {
	for(unsigned int i=0; i<N; i++) {
		if(x[i] < rhs.x[i])
			return true;
		if(x[i] > rhs.x[i])
			return false;
	}
	return false;
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator<=(const Vec<N,T>& rhs) const {
	for(unsigned int i=0; i<N; i++) {
		if(x[i] < rhs.x[i])
			return true;
		if(x[i] > rhs.x[i])
			return false;
	}
	return true;
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator>(const Vec<N,T>& rhs) const {
	for(unsigned int i=0; i<N; i++) {
		if(x[i] > rhs.x[i])
			return true;
		if(x[i] < rhs.x[i])
			return false;
	}
	return false;
}

template<unsigned int N, typename T>
bool Vec<N,T>::operator>=(const Vec<N,T>& rhs) const {
	for(unsigned int i=0; i<N; i++) {
		if(x[i] > rhs.x[i])
			return true;
		if(x[i] < rhs.x[i])
			return false;
	}
	return true;
}

/// string output representation for vectors
template<unsigned int N, typename T>
std::ostream& operator<<(std::ostream& o, const Vec<N,T>& v) {
	o << "<\t";
	for(unsigned int i=0; i<N; i++) {
		if(i) o << ",\t";
		o << v[i];
	}
	o << "\t>";
	return o;
}


template<unsigned int N, typename T>
void Vec<N,T>::display(const char* suffix) const {
	std::cout << *this << suffix;
}


/// cross product of 2-vectors
template<typename T>
T cross( const Vec<2,T>& v1, const Vec<2,T>& v2 ) { return v1[0]*v2[1]-v1[1]*v2[0]; }

/// cross product of 3-vectors
template<typename T>
Vec<3,T> cross( const Vec<3,T>& a, const Vec<3,T>& b ) {
	return Vec<3,T>(a[1]*b[2]-b[1]*a[2], a[2]*b[0]-b[2]*a[0], a[0]*b[1]-b[0]*a[1]);
}

/// rotation of a 2-vector 90 degrees counterclockwise
template<typename T>
Vec<2,T> rhOrtho( const Vec<2,T>& v ) {
	return Vec<2,T>(-v[1],v[0]); 
}

/// rotation of a 2-vector by given angle
template<typename T>
Vec<2,T> rotated(const Vec<2,T>& v, T a ) {
	return Vec<2,T>( v[0]*cos(a)-v[1]*sin(a), v[1]*cos(a)+v[0]*sin(a) );
}

/// orthonormal 2-vector 90 degrees counterclockwise of given 2-vector
template<typename T>
Vec<2,T> rhOrthoNorm( const Vec<2,T>& v ) { 
	return Vec<2,T>(-v[1],v[0]).normalized();
}

/// atan2() angle of a 2-vector
template<typename T>
T angle( const Vec<2,T>& v ) { 
	return atan2(v[1],v[0]);
}
 
/// vec2 from polar form specification
template<typename T>
Vec<2,T> polarVec(T r, T th) {
	return Vec<2,T>(r*cos(th),r*sin(th));
}

/// Vec to vector<double>
template<unsigned int N, typename T>
std::vector<double> vec2doublevec(const Vec<N,T>& v) {
	std::vector<double> dv(N);
	for(unsigned int i=0; i<N; i++) dv[i] = (double)v[i];
	return dv;
}


#endif
