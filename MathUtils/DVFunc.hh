#ifndef DVFUNC_HH
#define DVFUNC_HH 1

#include "Vec.hh"
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

/// virtual base class for a vector-valued differential function
template<unsigned int N1, unsigned int N2, typename T>
class DVFunc {
public:
	/// constructor
	DVFunc() {}
	/// destructor
	virtual ~DVFunc() {}
	
	/// evaluate function
	virtual Vec<N2,T> operator()(const Vec<N1,T>& x) const = 0;
	
	/// partial derivative along axis i
	virtual Vec<N2,T> deriv(const Vec<N1,T>& x, unsigned int i) const;
};

template<unsigned int N1, unsigned int N2, typename T>
Vec<N2,T> DVFunc<N1,N2,T>::deriv(const Vec<N1,T>& x, unsigned int i) const {
	double h = 1e-6;
	Vec<N1,T> dh;
	dh[i] += h;
	return ((*this)(x+dh)-(*this)(x-dh)) / (2.*h);
}

/// virtual base class for a vector-valued differential function of 1 variable
template<unsigned int N, typename T>
class DVFunc1 {
public:
	/// constructor
	DVFunc1(): period(0) {}
	/// destructor
	virtual ~DVFunc1() {}
	
	T period;	//< period of function, if periodic
		
	/// evaluate function
	virtual Vec<N,T> operator()(T x) const = 0;
	
	/// derivative
	virtual Vec<N,T> deriv(T x) const;
	
	/// verify derivative calculation
	void check_deriv() const {
		for(int i=-10; i<21; i++)
			std::cout << 0.1*i << ": " << deriv(0.1*i) << " " << DVFunc1<N,T>::deriv(0.1*i) << std::endl;
	}
};

template<unsigned int N, typename T>
Vec<N,T> DVFunc1<N,T>::deriv(T x) const {
	double h = 1e-6;
	return ((*this)(x+h)-(*this)(x-h)) / (2.*h);
}

/// virtual base class for a differentiable univariate function
template<typename T>
class DFunc {
public:
	/// constructor
	DFunc() {}
	/// destructor
	virtual ~DFunc() {}
	
	/// evaluate function
	virtual T operator()(T x) const = 0;
	
	/// partial derivative along axis i
	virtual T deriv(T x) const;
};

template<typename T>
T DFunc<T>::deriv(T x) const {
	double h = 1e-6;
	return ((*this)(x+h)-(*this)(x-h)) / (2.*h);
}


#endif
