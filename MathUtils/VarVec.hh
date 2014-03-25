#ifndef VARVEC_HH
#define VARVEC_HH 1

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <vector>
#include "Vec.hh"
#include "Permutation.hh"
#include <iostream>
#include <algorithm>
#include <exception>

class DimensionMismatchError: public std::exception {
	virtual const char* what() const throw() {
		return "Dimension mismatch error!";
	}
};

/// Dynamically allocated length vectors
template<class T>
class VarVec {
public:
	/// Constructor
	VarVec(unsigned int n = 0): data(n) {}
	/// Constructor with fill element
	VarVec(unsigned int n, const T& i): data(n,i) {}
	/// Constructor from fixed-length vector
	template<unsigned int N>
	VarVec(const Vec<N,T>& v): data(&v[0],&v[0]+N) { }
	/// construct from start/end pointers
	template <class InputIterator>
	VarVec(InputIterator first, InputIterator last): data(first,last) {}
	/// Destructor
	~VarVec() {}
	
	/// mutable element access operator
	T& operator[](unsigned int i) { assert(i<data.size()); return data[i]; }
	/// immutable element access operator
	const T& operator[](unsigned int i) const { assert(i<data.size()); return data[i]; }
	/// immutable access to the whole data vector
	const std::vector<T>& getData() const { return data; }
	/// mutable access to the whole data vector
	std::vector<T>& getData() { return data; }
	/// pointer to beginning of array
	T* getDataPtr() { return &data.front(); }
	/// pointer to beginning of array
	const T* getDataPtr() const { return &data.front(); }
	/// append
	void push_back(const T& x) { data.push_back(x); }
	/// generate sub-vector
	VarVec<T> subvec(unsigned int a, unsigned int b) const { VarVec<T> V; V.data = std::vector<T>(&data[a],&data[b]); return V; }
	/// copy data from a sub-vector, starting at position i
	void load_subvec(const VarVec<T>& V, unsigned int i) { assert(i+V.size()<=size()); std::copy(V.getData().begin(), V.getData().end(), &data[i]); }
	
	/// size of vector
	unsigned int size() const { return data.size(); }
	/// display to stdout
	void display() const { std::cout << "< "; for(unsigned int i=0; i<size(); i++) std::cout << data[i] << " "; std::cout << ">\n"; }
	
	/// throw error if dimensions mismatches
	void checkDimensions(const VarVec& v) const throw(DimensionMismatchError) { if(v.size() != size()) throw(DimensionMismatchError()); }
	
	/// dot product with another vector
	T dot(const VarVec<T>& v) const {
		checkDimensions(v);
		if(!size()) return 0;
		T s = data[0]*v[0];
		for(unsigned int i=1; i<size(); i++) s+=data[i]*v[i];
		return s;
	}
	/// square magnitude \f$ v \cdot v \f$
	T mag2() const { return dot(*this); }
	/// magnitude \f$ \sqrt{v\cdot v} \f$
	T mag() const { return sqrt(mag2()); }
	/// L2 norm of this vector (= mag() for float/double)
	double norm_L2() const;
	/// maximum L2-norm over all elements
	double max_norm_L2() const;
	/// sum of vector elements
	T sum() const { T s = data[0]; for(unsigned int i=1; i<size(); i++) s += data[i]; return s; }
	/// product of vector elements
	T prod() const { T s = data[0]; for(unsigned int i=1; i<size(); i++) s *= data[i]; return s; }
	/// minimum element of vector
	T min() const { return *std::min_element(data.begin(),data.end()); }
	/// maximum element of vector
	T max() const { return *std::max_element(data.begin(),data.end()); }
		
	/// this vector, normalized to magnitude 1
	VarVec<T> normalized() const { return (*this)/mag(); }
	/// project out component parallel to another vector
	VarVec<T> paraProj(const VarVec<T>& v) const { return v*(dot(v)/v.mag2()); }
	/// project out component orthogonal to another vector
	VarVec<T> orthoProj(const VarVec<T>& v) const { return (*this)-paraProj(v); }
	/// angle with another vector
	//T angle(const VarVec<T> v) const { return acos(dot(v)/sqrt(mag2()*v.mag2())); }
	
	
	/// unary minus operator
	const VarVec<T> operator-() const;
	
	/// inplace addition
	VarVec<T>& operator+=(const VarVec<T>& rhs);	
	/// inplace subtraction
	VarVec<T>& operator-=(const VarVec<T>& rhs);	
	/// inplace addition
	VarVec<T>& operator+=(T c);	
	/// inplace subtraction
	VarVec<T>& operator-=(T c);
	
	/// inplace multiplication
	VarVec<T>& operator*=(T c);
	/// inplace elementwise multiplication
	VarVec<T>& operator*=(const VarVec<T>& other);
	/// inplace division
	VarVec<T>& operator/=(T c);
	/// inplace elementwise division
	VarVec<T>& operator/=(const VarVec<T>& other);
	
	/// addition operator
	const VarVec<T> operator+(const VarVec<T>& other) const;	
	/// subtraction operator
	const VarVec<T> operator-(const VarVec<T>& other) const;	
	/// addition operator
	const VarVec<T> operator+(T c) const;	
	/// subtraction operator
	const VarVec<T> operator-(T c) const;
	
	/// multiplication operator
	const VarVec<T> operator*(T c) const;
	/// elementwise multiplication operator
	const VarVec<T> operator*(const VarVec<T>& other) const;
	/// division operator
	const VarVec<T> operator/(T c) const;
	/// elementwise division operator
	const VarVec<T> operator/(const VarVec<T>& other) const;
	
	/// equality operator
	bool operator==(const VarVec<T>& rhs) const;
	/// inequality operator
	bool operator!=(const VarVec<T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator<(const VarVec<T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator<=(const VarVec<T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator>(const VarVec<T>& rhs) const;
	/// dictionary-order comparison operator
	bool operator>=(const VarVec<T>& rhs) const;
	
	/// zero all elements
	VarVec& zero() { for(unsigned int i=0; i<size(); i++) data[i] = T(); return *this; }
	/// Make the nth element of the vector =1, all others =0
	VarVec<T>& basis(int n) { zero(); data[n] += 1.0; return *this; }
	/// Fill the vector with random numbers in [0,1]
	VarVec<T>& random() { zero(); for(unsigned int i=0; i<size(); i++) data[i] += (double)rand()/(double)RAND_MAX; return *this; }
	/// Fill the vector with ascending sequence \f$ r_0, r_0+1, r_0+2, \cdots \f$
	VarVec<T>& ramp(T r0) { for(unsigned int i=0; i<size(); i++) data[i] = r0 + i; return *this; }

	/// Create a new Vec by permuting the order of this vector's elements
	const VarVec<T> permuted(const Permutation& p) const;
	/// Permute the order of this vector's elements
	VarVec<T>& permute(const Permutation& p);
	
protected:
	std::vector<T> data;
	
};


template<typename T>
const VarVec<T> VarVec<T>::operator-() const {
	VarVec<T> v = *this;
	for(unsigned int i=0; i<size(); i++)
		v[i] *= -1.0;
	return v;
}

template<typename T>
VarVec<T>& VarVec<T>::operator+=(const VarVec<T>& rhs) {
	checkDimensions(rhs);
	for(unsigned int i=0; i<size(); i++)
		data[i] += rhs[i];
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator-=(const VarVec<T>& rhs) {
	checkDimensions(rhs);
	for(unsigned int i=0; i<size(); i++)
		data[i] -= rhs[i];
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator+=(T c) {
	for(unsigned int i=0; i<size(); i++)
		data[i] += c;
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator-=(T c) {
	for(unsigned int i=0; i<size(); i++)
		data[i] -= c;
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator*=(T c) {
	for(unsigned int i=0; i<size(); i++)
		data[i] *= c;
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator*=(const VarVec<T>& other) {
	checkDimensions(other);
	for(unsigned int i=0; i<size(); i++)
		data[i] *= other[i];
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator/=(T c) {
	for(unsigned int i=0; i<size(); i++)
		data[i] /= c;
	return *this;
}

template<typename T>
VarVec<T>& VarVec<T>::operator/=(const VarVec<T>& other) {
	checkDimensions(other);
	for(unsigned int i=0; i<size(); i++)
		data[i] /= other[i];
	return *this;
}

template<typename T>
const VarVec<T> VarVec<T>::operator+(const VarVec<T>& other) const {
	VarVec<T> result = *this;
	result += other;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator-(const VarVec<T>& other) const {
	VarVec<T> result = *this;
	result -= other;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator+(T c) const {
	VarVec<T> result = *this;
	result += c;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator-(T c) const {
	VarVec<T> result = *this;
	result -= c;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator*(T c) const {
	VarVec<T> result = *this;
	result *= c;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator*(const VarVec<T>& other) const {
	VarVec<T> result = *this;
	result *= other;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator/(T c) const {
	VarVec<T> result = *this;
	result /= c;
	return result;
}

template<typename T>
const VarVec<T> VarVec<T>::operator/(const VarVec<T>& other) const {
	VarVec<T> result = *this;
	result /= other;
	return result;
}

 
template<typename T>
bool VarVec<T>::operator==(const VarVec<T>& rhs) const {
	checkDimensions(rhs);
	for(unsigned int i=0; i<size(); i++)
		if(data[i] != rhs.data[i])
			return false;
	return true;
}

template<typename T>
bool VarVec<T>::operator!=(const VarVec<T>& rhs) const {
	return !(*this == rhs);
}

template<typename T>
bool VarVec<T>::operator<(const VarVec<T>& rhs) const {
	checkDimensions(rhs);
	for(unsigned int i=0; i<size(); i++) {
		if(data[i] < rhs.data[i])
			return true;
		if(data[i] > rhs.data[i])
			return false;
	}
	return false;
}

template<typename T>
bool VarVec<T>::operator<=(const VarVec<T>& rhs) const {
	checkDimensions(rhs);
	for(unsigned int i=0; i<size(); i++) {
		if(data[i] < rhs.data[i])
			return true;
		if(data[i] > rhs.data[i])
			return false;
	}
	return true;
}

template<typename T>
bool VarVec<T>::operator>(const VarVec<T>& rhs) const {
	return !((*this) <= rhs);
}

template<typename T>
bool VarVec<T>::operator>=(const VarVec<T>& rhs) const {
	return !((*this) < rhs);
}


template<class T>
const VarVec<T> VarVec<T>::permuted(const Permutation& p) const
{
	VarVec<T> o = VarVec<T>(size());
	for(unsigned int i=0; i<size(); i++) o[i] = data[p[i]];
	return o;
}

template<class T>
VarVec<T>& VarVec<T>::permute(const Permutation& p)
{
	std::vector<T> dnew = std::vector<T>(size());
	for(unsigned int i=0; i<size(); i++) dnew[i] = data[p[i]];
	data = dnew;
	return *this;
}

/// VarVec to vector<double>
template<typename T>
std::vector<double>
varvec2doublevec(const VarVec<T>& v) {
	std::vector<double> dv(v.size());
	for(unsigned int i=0; i<v.size(); i++) dv[i] = (double)v[i];
	return dv;
}

namespace VarVec_element_norm_L2 {
	template<typename T>
	inline double norm_L2(const T& t) { return t.norm_L2(); }
  	template<>
	inline double norm_L2(const float& t) { return fabs(t); }
	template<>
	inline double norm_L2(const double& t) { return fabs(t); }
}

template<typename T>
double VarVec<T>::max_norm_L2() const {
	std::vector<double> vn;
	for(typename std::vector<T>::const_iterator it = data.begin(); it != data.end(); it++)
		vn.push_back(VarVec_element_norm_L2::norm_L2(*it));
	return *std::max_element(vn.begin(), vn.end());
}

template<typename T>
double VarVec<T>::norm_L2() const {
	double s = 0;
	for(typename std::vector<T>::const_iterator it = data.begin(); it != data.end(); it++) {
		double n = VarVec_element_norm_L2::norm_L2(*it);
		s += n*n;
	}
	return sqrt(s);
}

/// string output representation for vectors
template<typename T>
std::ostream& operator<<(std::ostream& o, const VarVec<T>& v) {
	o << "< ";
	for(unsigned int i=0; i<v.size(); i++)
		o << v[i] << " ";
	o << ">";
	return o;
}

#endif
