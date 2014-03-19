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
	DVFunc1() {}
	/// destructor
	virtual ~DVFunc1() {}

	/// evaluate function
	virtual Vec<N,T> operator()(T x) const = 0;
	
	/// derivative
	virtual Vec<N,T> deriv(T x) const;
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



//
//==========================================================================
//


/// append sequential univariate functions, assumed to each reside on [0,1]
template<unsigned int N, typename T>
class PathJoiner: public DVFunc1<N,T> {
public:
	/// constructor
	PathJoiner() {
		seg_endpts.push_back(0);
	}
	
	/// destructor
	virtual ~PathJoiner() {}
	
	/// evaluate function
	virtual Vec<N,T> operator()(T x) const {
		unsigned int i = locate(x);
		return seg_offsets[i] + (*mysegs[i])(x);
	}
	
	/// derivative
	virtual Vec<N,T> deriv(T x) const {
		unsigned int i = locate(x);
		return mysegs[i]->deriv(x) * seg_endpts.back()/(seg_endpts[i+1]-seg_endpts[i]);
	}
	
	/// append a segment, spaced to maintain continuous pathlength derivative
	virtual void append( DVFunc1<N,T>* f, bool smooth_pathlength = true) {
		assert(f);
		if(!mysegs.size()) {
			seg_endpts.push_back(1);
		} else {
			unsigned int i = mysegs.size();
			T mul = smooth_pathlength ? f->deriv(0).mag() / mysegs.back()->deriv(1).mag() * (seg_endpts[i]-seg_endpts[i-1]) : 1;
			seg_endpts.push_back(seg_endpts.back() + mul);
		}
		if(mysegs.size())
			seg_offsets.push_back(seg_offsets.back() + (*mysegs.back())(1.) - (*f)(0));
		else
			seg_offsets.push_back(Vec<N,T>());
		mysegs.push_back(f);
	}
	
	/// display contents
	void display() const {
		std::cout << "PathJoiner in " << mysegs.size() << " segments:\n";
		for(unsigned int i=0; i<mysegs.size(); i++)
			std::cout << "\t" << seg_endpts[i] << " to " << seg_endpts[i+1] << "\tstarting at " << seg_offsets[i] << std::endl;
	}
	
protected:

	/// locate sub-function and position therein
	unsigned int locate(T& l) const {
		assert(mysegs.size());
		typename std::vector<T>::const_iterator it = std::upper_bound(seg_endpts.begin(), seg_endpts.end(), l*seg_endpts.back());
		unsigned int i = it-seg_endpts.begin();
		if(!i) i += 1;
		if(it==seg_endpts.end()) i -= 1;
		l = (l*seg_endpts.back()-seg_endpts[i-1])/(seg_endpts[i]-seg_endpts[i-1]);
		return i-1;
	}

	/// delete all items
	void clear() {
		for(unsigned int i=0; i<mysegs.size(); i++)
			delete mysegs[i];
		mysegs.clear();
		seg_endpts.clear();
		seg_offsets.clear();
		seg_endpts.push_back(0);
	}
	
	std::vector< DVFunc1<N,T>* > mysegs;	//< subsections
	std::vector<T> seg_endpts;				//< map from global l to sub-range endpoints
	std::vector< Vec<N,T> > seg_offsets;	// endpoint offsets for each segment
};



#endif
