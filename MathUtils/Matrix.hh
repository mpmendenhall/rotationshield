#ifndef MATRIX_HH
#define MATRIX_HH 1

#include "Vec.hh"
#include <vector>

/// A templatized, fixed size, statically allocated matrix class.
/**
 Not particularly optimized or clever, but convenient for smallish matrices
 or matrices of unusual special types (e.g. a matrix of circulant matrices)
 template<unsigned int N, unsigned int M, typename T>
*/
template<unsigned int N, unsigned int M, typename T>
class Matrix {
public:
	/// constructor
	Matrix(): vv(Vec<N*M,T>()) { }
	
	/// constructor from vector
	Matrix(const Vec<N*M,T>& v): vv(v) {}
	
	/// destructor
	~Matrix() {}

	/// generate a random-filled matrix
	static Matrix<N,M,T> random();
	
	/// const element access
	const T& operator()(unsigned int n, unsigned int m) const { assert(n<N && m<M); return vv[n*M+m]; }
	/// mutable element access
	T& operator()(unsigned int n, unsigned int m) { assert(n<N && m<M); return vv[n*M+m]; }
	/// direct access to data vector
	const Vec<N*M,T>& getData() const { return vv; }
	/// mutable vector element access
	T& operator[](unsigned int i) { return vv[i]; }
	/// const vector element access
	const T& operator[](unsigned int i) const { return vv[i]; }
	/// row vector
	Vec<M,T> row(unsigned int i) const;
	/// column vector
	Vec<N,T> col(unsigned int i) const;
	
	/// transpose
	Matrix<M,N,T> transposed() const;
	///unary minus
	const Matrix<N,M,T> operator-() const;
	
	/// inplace multiplication by a constant
	void operator*=(const T& c) { vv *= c; }
	/// multiplication by a constant
	const Matrix<N,M,T> operator*(const T& c) const;
	/// inplace division by a constant
	void operator/=(const T& c) { vv /= c; }
	/// division by a constant
	const Matrix<N,M,T> operator/(const T& c) const;
	/// inplace addition of a matrix
	void operator+=(const Matrix<N,M,T>& rhs) { vv += rhs.getData(); }
	/// addition of a matrix
	const Matrix<N,M,T> operator+(const Matrix<N,M,T>& rhs) const;
	/// inplace subtraction of a matrix
	void operator-=(const Matrix<N,M,T>& rhs) { vv -= rhs.getData(); }
	/// subtraction of a matrix
	const Matrix<N,M,T> operator-(const Matrix<N,M,T>& rhs) const;
	
	/// matrix multiplication
	template<unsigned int L>
	const Matrix<N,L,T> operator*(const Matrix<M,L,T>& B) const;
	/// left-multiply a vector
	const Vec<N,T> lMultiply(const Vec<M,T>& v) const;
	/// right-multiply a vector
	const Vec<M,T> rMultiply(const Vec<N,T>& v) const;
	/// vector multiplication
	const Vec<N,T> operator*(const Vec<M,T>& v) const { return lMultiply(v); }
	
	/// matrix multiplication and assignment
	void operator*=(const Matrix<M,M,T>& B) { (*this) = B*(*this); }
	
	/// inplace inverse
	const Matrix<N,M,T>& invert();
	
private:
	Vec<N*M,T> vv;
	
	/// step in inversion process
	void subinvert(unsigned int n);
};

template<unsigned int N, unsigned int M, typename T>
Matrix<N,M,T> Matrix<N,M,T>::random() { 
	Matrix<N,M,T> foo; 
	for(unsigned int i=0; i<N*M; i++)
		foo[i] = 0.1+T(rand())/T(RAND_MAX);
	return foo; 
}


template<unsigned int N, unsigned int M, typename T>
Matrix<M,N,T> Matrix<N,M,T>::transposed() const {
	Matrix<M,N,T> foo;
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			foo(c,r) = (*this)(r,c);
	return foo;
}

template<unsigned int N, unsigned int M, typename T>
Vec<M,T> Matrix<N,M,T>::row(unsigned int i) const {
	Vec<M,T> v;
	for(unsigned int j=0; j<M; j++) v[i] = (*this)(i,j);
	return v;
}

template<unsigned int N, unsigned int M, typename T>
Vec<N,T> Matrix<N,M,T>::col(unsigned int i) const {
	Vec<N,T> v;
	for(unsigned int j=0; j<N; j++) v[i] = (*this)(j,i);
	return v;
}


template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T> Matrix<N,M,T>::operator-() const {
	Matrix<N,M,T> foo; 
	for(unsigned int i=0; i<N*M; i++)
		foo[i] = -(*this)[i];
	return foo;
}

template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T> Matrix<N,M,T>::operator*(const T& c) const {
	Matrix<N,M,T> foo = *this;
	foo *= c;
	return foo;
}

template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T> Matrix<N,M,T>::operator/(const T& c) const {
	Matrix<N,M,T> foo = *this;
	foo /= c;
	return foo;
}

template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T> Matrix<N,M,T>::operator+(const Matrix<N,M,T>& rhs) const {
	Matrix<N,M,T> foo = *this;
	foo += rhs;
	return foo;
}


template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T> Matrix<N,M,T>::operator-(const Matrix<N,M,T>& rhs) const {
	Matrix<N,M,T> foo = *this;
	foo -= rhs;
	return foo;
}

template<unsigned int N, unsigned int M, typename T>
template<unsigned int L>
const Matrix<N,L,T> Matrix<N,M,T>::operator*(const Matrix<M,L,T>& B) const {
	Matrix<N,L,T> C = Matrix<N,L,T>();
	for(unsigned int r=0; r<N; r++) {
		for(unsigned int c=0; c<L; c++) {
			C(r,c) = (*this)(r,0)*B(0,c);
			for(unsigned int i=1; i<M; i++)
				C(r,c) += (*this)(r,i)*B(i,c);
		}
	}
	return C;
}

template<unsigned int N, unsigned int M, typename T>
const Vec<N,T> Matrix<N,M,T>::lMultiply(const Vec<M,T>& v) const {
	Vec<N,T> a = Vec<N,T>();
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			a[r] += (*this)(r,c)*v[c];
	return a;
}

template<unsigned int N, unsigned int M, typename T>
const Vec<M,T> Matrix<N,M,T>::rMultiply(const Vec<N,T>& v) const {
	Vec<M,T> a = Vec<M,T>();
	for(unsigned int r=0; r<N; r++)
		for(unsigned int c=0; c<M; c++)
			a[c] += v[r]*(*this)(r,c);
	return a;
}


template<unsigned int N, unsigned int M, typename T>
const Matrix<N,M,T>& Matrix<N,M,T>::invert() {
	assert(N==M);
	subinvert(0);
	return *this;
}





namespace matrix_element_inversion {
	
	template<typename T>
	inline void invert_element(T& t) { t.invert(); }
	template<>
	inline void invert_element(float& t) { t = 1.0/t; }
	template<>
	inline void invert_element(double& t) { t = 1.0/t; }
	
}

template<unsigned int N, unsigned int M, typename T>
void Matrix<N,M,T>::subinvert(unsigned int n) {
	
	// invert the first cell
	T& firstcell = (*this)(n,n);
	matrix_element_inversion::invert_element(firstcell);
	for(unsigned int i=n+1; i<N; i++)
		(*this)(n,i) *= firstcell;
		//(*this)(n,i) = firstcell*(*this)(n,i);
	
	// use to clear first column
	for(unsigned int r=n+1; r<N; r++) {
		T& m0 = (*this)(r,n);
		for(unsigned int c=n+1; c<N; c++)
			(*this)(r,c) -= (*this)(n,c)*m0;
		m0 *= -firstcell;
		//m0 = -m0*firstcell;
	}
	
	if(n==N-1)
		return;
	
	//invert the submatrix
	subinvert(n+1);
	
	// temporary space allocation
	std::vector<T> subvec = std::vector<T>(N-n-1);
	
	// first column gets multiplied by inverting submatrix
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



/// string output representation for matrix
template<unsigned int N, unsigned int M, typename T>
std::ostream& operator<<(std::ostream& o, const Matrix<N,M,T>& A) {
	for(unsigned int r=0; r<N; r++) {
		o << "| ";
		for(unsigned int c=0; c<M; c++)
			o << A(r,c) << " ";
		o << "|\n";
	}
	return o;
}

#endif
