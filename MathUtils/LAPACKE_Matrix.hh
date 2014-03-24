#ifndef LAPACKE_MATRIX_HH
#define LAPACKE_MATRIX_HH 1

#include <complex.h>
#define lapack_complex_float    float _Complex
#define lapack_complex_double   double _Complex
#include <lapacke.h>

#include "ComplexT.hh"
#include <algorithm>
#include <cassert>

typedef ComplexT<double>  complex<double> ;
typedef int MKL_INT;

/// Templatized wrapper for LAPACKE matrix operations; in column-major form
template<typename T>
class LAPACKE_Matrix {
public:
	/// constructor
	LAPACKE_Matrix(unsigned int m, unsigned int n): M(m), N(n) { data = new T[N*M]; }
	/// destructor
	~LAPACKE_Matrix() { delete[] data; }
	
	/// const element access
	T operator()(unsigned int n, unsigned int m) const { assert(n<N && m<M); return data[n+m*N]; }
	/// mutable element access
	T& operator()(unsigned int n, unsigned int m) { assert(n<N && m<M); return data[n+m*N]; }
	
	const unsigned int M;	//< number of rows
	const unsigned int N;	//< number of columns
	T* data;	//< data, in column-major order by default
		
	/// appropriate data type version of cblas multiply routine gemm, C = alpha*op(A)*op(B) + beta*C
	static void (*f_gemm)(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE,
						  const MKL_INT, const MKL_INT, const MKL_INT, const double, const T*,
						  const MKL_INT, const T*, const MKL_INT, double, T*, const MKL_INT);

};

//void (*)(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, const MKL_INT, const MKL_INT, const MKL_INT, const void *, const void *, const MKL_INT, const void *, const MKL_INT, const void *, void *, const MKL_INT)
//void (*)(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, const int,     const int,     const int,     const double, const double *, const int, const double *, const int, const double, double *, const int)': type mismatch at 7th parameter


template<>
void (*LAPACKE_Matrix<double>::f_gemm)(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE,
									   const MKL_INT, const MKL_INT, const MKL_INT, const double, const double*,
									   const MKL_INT, const double*, const MKL_INT, double, double*, const MKL_INT) = &cblas_dgemm;
template<>
void (*LAPACKE_Matrix< complex<double> >::f_gemm)(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE,
										const MKL_INT, const MKL_INT, const MKL_INT, const double, const  complex<double> *,
										const MKL_INT, const  complex<double> *, const MKL_INT, double,  complex<double> *, const MKL_INT) = &cblas_zgemm;

/// form a new matrix by multiplication C = A*B
template<typename T>
LAPACKE_Matrix* operator*(const LAPACKE_Matrix<T>& A, const LAPACKE_Matrix<T>& B) {
	assert(A.N == B.M);
	
	LAPACKE_Matrix<T>* C = new LAPACKE_Matrix<T>(A.M,B.N);
	
	const T alpha(1.);
	const T beta(0.);
	(*LAPACKE_Matrix<T>::f_gemm)(CblasColMajor,	// data order, CblasColMajor or CblasRowMajor
								 'n',			// char for op(A), 'n'=none, 't'=transpose, 'c'=conjugate
								 'n',			// char for op(B), 'n'=none, 't'=transpose, 'c'=conjugate
								 C.M,			// rows of op(A) and C
								 C.N,			// columns of op(B) and C
								 A.N,			// columns of op(A), rows of op(B)
								 alpha,			// pointer to scalar alpha
								 A.data,		// matrix A's data
								 A.M,			// leading dimension of op(A)
								 B.data,		// matrix B's data
								 B.M,			// B's leading dimension
								 beta,			// pointer to scalar beta
								 C.data,		// output matrix C data
								 C.M			// leading dimension of C
								 );
	return C;
}




/// Templatized wrapper for SVD of matrix A = U S V^T
template<typename T>
class LAPACKE_Matrix_SVD {
	/// Constructor
	LAPACKE_Matrix_SVD(const LAPACKE_Matrix<T>& A): S(std::min(A.M,A.N),1), U(A.M,S.M), VT(S.M,A.N) {
		
		lapack_int info;
		
		char diag = (A.M >= A.N ? 'U':'L');		// upper or lower diagonal reduction, depending on A's dimensions
		unsigned int dimB = std::min(A.M,A.N);	// dimension of bidiagonal matrix B
		double* e = new double[S.M];
		T* tauQ = new T[S.M];
		T* tauP = new T[S.M];
		
		// overall reduce general A = (U*Q) * S * (P^H * V^T)
		
		// decomposes m*n A = Q B P^H,  Q and P unitary, B bidiagonal of min(m,n)*min(m,n)
		// over-writes A with B on diagonals, P above, and Q below.
		info = (*f_gebrd)(LAPACK_COL_MAJOR,	// LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR data ordering,
						  A.M,				// m: number of rows in A
						  A.N,				// n: number of columns in A
						  A.data,			// matrix A
						  A.M,				// lda: leading dimension of A, >= max(1,m)
						  S.data,			// diagonal elements of B, dimension >= max(1, min(m, n))
						  e,				// off-diagonal elements of B, dimension >= max(1, min(m, n) - 1)
						  tauQ,				// array with extra into on Q, dimension >= max(1, min(m, n))
						  tauP				// array with extra info on P, dimension >= max(1, min(m, n))
						  );
		assert(!info);
		delete[] tauQ;
		delete[] tauP;
		
		// TODO transfer data into new matrices!
		LAPACKE_Matrix<T> U();
		LAPACKE_Matrix<T> VT();
		for(unsigned int m=0; m<A.M; m++) {
			for(unsigned int n=0; n<A.N; n++) {
			}
		}
		assert(false);
		
		// decompose B = U2 S V2^H, Q and P orthogonal singular vectors, S diagonal singular values
		// overwrites d = S; destroys e
		info = (*f_bdsqr)(LAPACK_COL_MAJOR,	// LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR data ordering
						  diag,				// 'U'pper or 'L'ower bidiagonal
						  S.M,				// order >= 0 of matrix B
						  V.N,				// ncvt: number of columns of V^T, right singular vectors to calculate
						  U.M,				// nru: number of rows of U, left singular vectors to calculate
						  0,				// ncc: number of columns in C used for QH*C; set 0 if no "C" supplied
						  S.data,			// d: diagonal elements of B; must be size >= max(1,n)
						  e,				// e: n-1 off-diagonal elements of B; must be size >= max(1,n)
						  VT.data,			// n * ncvt matrix for right singular vectors; unused if ncvt=0
						  VT.M,				// leading dimension of VT; >= max(1,n) for ncvt>0; >= 1 otherwise.
						  U.data,			// nru * n unit matrix U; unused if nru=0
						  U.M,				// leading dimension of U; >= max(1,nru)
						  NULL,				// matrix for calculating Q^H*C; second dimension >= max(1,ncc); unused if ncc = 0
						  1					// leading dimension of C; >= max(1,n) if ncc>0; >=1 otherwise
						  );
		assert(!info);
		
		delete[] e;
	}
	
	/// solve linear least squares problem, minimizing || B - A*X ||_2 for X given B
	void lstsq_solve(LAPACKE_Matrix<T>& B) {
		lapack_int info;
		info = (*f_gelss)(LAPACK_COL_MAJOR,	// LAPACK_COL_MAJOR or LAPACK_ROW_MAJOR data ordering
						  U.M,				// rows of A
						  V.N,				// columns of A
						  B.N,				// number of right-hand-side columns in B
						  lapack_complex_double* a,
						  lapack_int lda,
						  lapack_complex_double* b,
						  lapack_int ldb,
						  double* s,
						  double rcond,
						  lapack_int* rank
						  );
	}
	
	LAPACKE_Matrix<T> S;		//< singular values diagonal
	LAPACKE_Matrix<T> U;		//< left singular vectors
	LAPACKE_Matrix<T> VT;		//< right
	
	/// appropriate data type version of xGEBRD bidiagonal reduction
	static lapack_int (*f_gebrd)(int, lapack_int, lapack_int, T*, lapack_int, double*, double*, T*, T*);
	/// appropriate data type version of zBDSQR SVD from bidiagonal
	static lapack_int (*f_bdsqr)(int, char, lapack_int, lapack_int, lapack_int, double*, double*, T*, lapack_int, T*, lapack_int, T*, lapack_int);
	/// appropriate data type Least Squares solution using SVD
	static lapack_int (*f_gelss)(int, lapack_int, lapack_int, lapack_int, T*, lapack_int, T*, lapack_int, double*, double, lapack_int*);
	
};

template<>
lapack_int (*LAPACKE_Matrix_SVD<double>::f_gebrd)(int, lapack_int, lapack_int, double*, lapack_int, double*, double*, double*, double*) = &LAPACKE_dgebrd;
template<>
lapack_int LAPACKE_Matrix_SVD< complex<double> >::f_gebrd = &LAPACKE_zgebrd;

template<>
lapack_int LAPACKE_Matrix_SVD<double>::f_bdsqr = &LAPACKE_dbdsqr;
template<>
lapack_int LAPACKE_Matrix_SVD< complex<double> >::f_bdsqr = &LAPACKE_zbdsqr;

template<>
lapack_int LAPACKE_Matrix_SVD<double>::f_gelss = &LAPACKE_dgelss;
template<>
lapack_int LAPACKE_Matrix_SVD< complex<double> >::f_gelss = &LAPACKE_zgelss;

#endif
