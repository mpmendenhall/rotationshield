#ifndef VARMAT_HH
#define VARMAT_HH 1

#include "VarVec.hh"
#include "Matrix.hh"
#include <vector>

/// A templatized, dynamically allocated matrix class.
/**
 Not particularly optimized or clever, but convenient for smallish matrices
 or matrices of unusual special types
 */
template<typename T>
class VarMat {
public:
	/// constructor
	VarMat(unsigned int n = 0, unsigned int m = 0): N(n), M(m), vv(VarVec<T>(n*m)) { }
	/// constructor from fixed matrix
	template<unsigned int NN, unsigned int MM>
	VarMat(Matrix<NN,MM,T> A): N(NN), M(MM), vv(A.getData()) {}
	/// destructor
	~VarMat() {}
	
	/// generate a random-filled VarMat
	static VarMat<T> random();
	
	/// get n rows
	unsigned int nRows() const { return N; }
	/// get m cols
	unsigned int nCols() const { return M; }
	/// get total size
	unsigned int size() const { return vv.size(); }
	
	/// const element access
	const T& operator()(unsigned int n, unsigned int m) const { assert(n<N && m<M); return vv[n*M+m]; }
	/// mutable element access
	T& operator()(unsigned int n, unsigned int m) { assert(n<N && m<M); return vv[n*M+m]; }
	/// direct access to data vector
	const VarVec<T>& getData() const { return vv; }
	/// mutable vector element access
	T& operator[](unsigned int i) { return vv[i]; }
	/// const vector element access
	const T& operator[](unsigned int i) const { return vv[i]; }
	
	/// transpose
	VarMat<T> transposed() const;
	///unary minus
	const VarMat<T> operator-() const;
	
	/// throw error if dimensions mismatches
	void checkDimensions(const VarMat& m) const throw(DimensionMismatchError) { if(m.nRows() != N || m.nCols() != M) throw(DimensionMismatchError()); }
	
	/// inplace multiplication by a constant
	void operator*=(const T& c) { vv *= c; }
	/// multiplication by a constant
	const VarMat<T> operator*(const T& c) const;
	/// inplace division by a constant
	void operator/=(const T& c) { vv /= c; }
	/// division by a constant
	const VarMat<T> operator/(const T& c) const;
	/// inplace addition of a VarMat
	void operator+=(const VarMat<T>& rhs) { checkDimensions(rhs); vv += rhs.getData(); }
	/// addition of a VarMat
	const VarMat<T> operator+(const VarMat<T>& rhs) const;
	/// inplace subtraction of a VarMat
	void operator-=(const VarMat<T>& rhs) { checkDimensions(rhs); vv -= rhs.getData(); }
	/// subtraction of a VarMat
	const VarMat<T> operator-(const VarMat<T>& rhs) const;
	
	/// VarMat multiplication
	const VarMat<T> operator*(const VarMat<T>& B) const;
	/// left-multiply a vector
	const VarVec<T> lMultiply(const VarVec<T>& v) const;
	/// right-multiply a vector
	const VarVec<T> rMultiply(const VarVec<T>& v) const;
	/// vector multiplication
	const VarVec<T> operator*(const VarVec<T>& v) const { return lMultiply(v); }
	
	/// VarMat multiplication and assignment
	void operator*=(const VarMat<T>& B) { (*this) = (*this)*B; }
	
	/// inplace inverse
	const VarMat<T>& invert();
	
private:
	
	unsigned int N;
	unsigned int M;
	VarVec<T> vv;
	
	/// step in inversion process
	void subinvert(unsigned int n);
};

template<typename T>
VarMat<T> VarMat<T>::random() { 
	VarMat<T> foo; 
	for(unsigned int i=0; i<foo.size(); i++)
		foo[i] = 0.1+T(rand())/T(RAND_MAX);
	return foo; 
}


template<typename T>
VarMat<T> VarMat<T>::transposed() const {
	VarMat<T> foo = VarMat(M,N);
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			foo(c,r) = (*this)(r,c);
	return foo;
}

template<typename T>
const VarMat<T> VarMat<T>::operator-() const {
	VarMat<T> foo = VarMat(N,M); 
	for(unsigned int i=0; i<N*M; i++)
		foo[i] = -(*this)[i];
	return foo;
}

template<typename T>
const VarMat<T> VarMat<T>::operator*(const T& c) const {
	VarMat<T> foo = *this;
	foo *= c;
	return foo;
}

template<typename T>
const VarMat<T> VarMat<T>::operator/(const T& c) const {
	VarMat<T> foo = *this;
	foo /= c;
	return foo;
}

template<typename T>
const VarMat<T> VarMat<T>::operator+(const VarMat<T>& rhs) const {
	VarMat<T> foo = *this;
	foo += rhs;
	return foo;
}


template<typename T>
const VarMat<T> VarMat<T>::operator-(const VarMat<T>& rhs) const {
	VarMat<T> foo = *this;
	foo -= rhs;
	return foo;
}

template<typename T>
const VarMat<T> VarMat<T>::operator*(const VarMat<T>& B) const {
	if(B.nRows() != M)
		throw(DimensionMismatchError());
	unsigned int L = B.nCols();
	VarMat<T> C = VarMat<T>(N,L);
	for(unsigned int r=0; r<N; r++) {
		for(unsigned int c=0; c<L; c++) {
			C(r,c) = (*this)(r,0)*B(0,c);
			for(unsigned int i=1; i<M; i++)
				C(r,c) += (*this)(r,i)*B(i,c);
		}
	}
	return C;
}

template<typename T>
const VarVec<T> VarMat<T>::lMultiply(const VarVec<T>& v) const {
	if(v.size() != M)
		throw(DimensionMismatchError());
	VarVec<T> a = VarVec<T>(N);
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			a[r] += (*this)(r,c)*v[c];
	return a;
}

template<typename T>
const VarVec<T> VarMat<T>::rMultiply(const VarVec<T>& v) const {
	if(v.size() != N)
		throw(DimensionMismatchError());
	VarVec<T> a = VarVec<T>(M);
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			a[c] += v[r]*(*this)(r,c);
	return a;
}

template<typename T>
const VarMat<T>& VarMat<T>::invert() {
	if(N != M)
		throw(DimensionMismatchError());
	subinvert(0);
	return *this;
}

namespace VarMat_element_inversion {
	
	template<typename T>
	inline void invert_element(T& t) { t.invert(); }
	template<>
	inline void invert_element(float& t) { t = 1.0/t; }
	template<>
	inline void invert_element(double& t) { t = 1.0/t; }
	
}

template<typename T>
void VarMat<T>::subinvert(unsigned int n) {
	
	// invert the first cell
	T& firstcell = (*this)(n,n);
	VarMat_element_inversion::invert_element(firstcell);
	for(unsigned int i=n+1; i<N; i++)
		(*this)(n,i) *= firstcell;
	
	// use to clear first column
	for(unsigned int r=n+1; r<N; r++) {
		T& m0 = (*this)(r,n);
		for(unsigned int c=n+1; c<N; c++)
			(*this)(r,c) -= (*this)(n,c)*m0;
		m0 *= -firstcell;
	}
	
	if(n==N-1)
		return;
	
	//invert the subVarMat
	subinvert(n+1);
	
	// temporary space allocation
	std::vector<T> subvec = std::vector<T>(N-n-1);
	
	// first column gets multiplied by inverting subVarMat
	for(unsigned int r=n+1; r<N; r++)
		subvec[r-n-1] = (*this)(r,n);
	for(unsigned int r=n+1; r<N; r++) {
		(*this)(r,n) = (*this)(r,n+1)*subvec[0];
		for(unsigned int c=n+2; c<N; c++)
			(*this)(r,n) += (*this)(r,c)*subvec[c-n-1];
	}
	
	//finish off by cleaning first row
	for(unsigned int c=n+1; c<N; c++)
		subvec[c-n-1] = (*this)(n,c);
	for(unsigned int c=n; c<N; c++) {
		if(c>n)
			(*this)(n,c) = -(*this)(n+1,c) * subvec[0];
		else
			(*this)(n,c) -= (*this)(n+1,c) * subvec[0];
		for(unsigned int r=n+2; r<N; r++)
			(*this)(n,c) -= (*this)(r,c) * subvec[r-n-1];
	}
	
}



/// string output representation for VarMat
template<typename T>
std::ostream& operator<<(std::ostream& o, const VarMat<T>& A) {
	for(unsigned int r=0; r<A.nRows(); r++) {
		o << "| ";
		for(unsigned int c=0; c<A.nCols(); c++)
			o << A(r,c) << " ";
		o << "|\n";
	}
	return o;
}

#endif
