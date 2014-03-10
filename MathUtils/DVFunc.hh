#ifndef DVFUNC_HH
#define DVFUNC_HH 1

#include "Vec.hh"

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

#endif
