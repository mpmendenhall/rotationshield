#ifndef COMPLEXT_HH
#define COMPLEXT_HH 1

#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <iostream>

template<typename T>
class ComplexT {
public:
	/// constructor (defaults to 0)
	ComplexT(T r = T(), T i = T()) { z[0] = r; z[1] = i; }
	
	/// square magnitude
	T mag2() const { return z[0]*z[0] + z[1]*z[1]; }
	/// magnitude
	T mag() const { return sqrt(mag2()); }

	/// unary minus operator
	const ComplexT<T> operator-() const;
	/// complex conjugate
	const ComplexT<T> conj() const { return ComplexT(z[0],-z[1]); }
	/// inverse 1/z
	const ComplexT<T> inverse() const;
	
	/// inplace addition
	ComplexT<T>& operator+=(const ComplexT<T>& rhs);	
	/// inplace subtraction
	ComplexT<T>& operator-=(const ComplexT<T>& rhs);	
	
	/// inplace multiplication
	ComplexT<T>& operator*=(T c);
	/// inplace elementwise multiplication
	ComplexT<T>& operator*=(const ComplexT<T>& other);
	/// inplace division
	ComplexT<T>& operator/=(T c);
	/// inplace elementwise division
	ComplexT<T>& operator/=(const ComplexT<T>& other);
	
	/// addition operator
	const ComplexT<T> operator+(const ComplexT<T>& other) const;	
	/// subtraction operator
	const ComplexT<T> operator-(const ComplexT<T>& other) const;	
	
	/// multiplication operator
	const ComplexT<T> operator*(T c) const;
	/// elementwise multiplication operator
	const ComplexT<T> operator*(const ComplexT<T>& other) const;
	/// division operator
	const ComplexT<T> operator/(T c) const;
	/// elementwise division operator
	const ComplexT<T> operator/(const ComplexT<T>& other) const;
	
	/// equality operator
	bool operator==(const ComplexT<T>& rhs) const;
	/// inequality operator
	bool operator!=(const ComplexT<T>& rhs) const;
	/// comparison operator
	bool operator<(const ComplexT<T>& rhs) const { return false; } //TODO
	/// comparison operator
	bool operator>(const ComplexT<T>& rhs) const { return false; } //TODO
	
	T z[2];
};


template<typename T>
const ComplexT<T> ComplexT<T>::operator-() const {
	return ComplexT(z[0],-z[1]);
}

template<typename T>
const ComplexT<T> ComplexT<T>::inverse() const {
	return conj()/mag2();
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator+=(const ComplexT<T>& rhs) {
	z[0] += rhs.z[0];
	z[1] += rhs.z[1];
	return *this;
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator-=(const ComplexT<T>& rhs) {
	z[0] -= rhs.z[0];
	z[1] -= rhs.z[1];
	return *this;
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator*=(T c) {
	z[0] *= c;
	z[1] *= c;
	return *this;
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator*=(const ComplexT<T>& other) {
	*this = ComplexT(z[0]*other.z[0] - z[1]*other.z[1], z[1]*other.z[0] + z[0]*other.z[1]);
	return *this;
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator/=(T c) {
	z[0] /= c;
	z[1] /= c;
	return *this;
}

template<typename T>
ComplexT<T>& ComplexT<T>::operator/=(const ComplexT<T>& other) {
	*this *= other.conj()/other.mag2();
	return *this;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator+(const ComplexT<T>& other) const {
	ComplexT<T> result = *this;
	result += other;
	return result;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator-(const ComplexT<T>& other) const {
	ComplexT<T> result = *this;
	result -= other;
	return result;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator*(T c) const {
	ComplexT<T> result = *this;
	result *= c;
	return result;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator*(const ComplexT<T>& other) const {
	ComplexT<T> result = *this;
	result *= other;
	return result;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator/(T c) const {
	ComplexT<T> result = *this;
	result /= c;
	return result;
}

template<typename T>
const ComplexT<T> ComplexT<T>::operator/(const ComplexT<T>& other) const {
	ComplexT<T> result = *this;
	result /= other;
	return result;
}

template<typename T>
bool ComplexT<T>::operator==(const ComplexT<T>& rhs) const {
	return z[0] == rhs.z[0] && z[1] == rhs.z[1];
}

template<typename T>
bool ComplexT<T>::operator!=(const ComplexT<T>& rhs) const {
	return !(*this == rhs);
}

/// string output representation for complex numbers
template<typename T>
std::ostream& operator<<(std::ostream& o, const ComplexT<T>& c) {
	o << c.z[0] << std::showpos << c.z[1] << "i";
	return o;
}

#endif
