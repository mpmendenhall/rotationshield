/// \file "CMatrix.hh" \brief Circulant matrices
#ifndef CMATRIX_HH
/// Make sure this header is only loaded once
#define CMATRIX_HH 1

#include <iostream>
#include <fftw3.h>
#include <vector>
#include <complex.h>
#include "VarVec.hh"

using namespace std;

/// Stores fftw data for FFT'ing 
struct cmatrix_fft {
	unsigned int ncyc; //< size of CMatrix this is meant for
	fftw_plan forwardplan; //< FFTW data for forward Fourier Transforms of this size
	fftw_plan reverseplan; //< FFTW data for inverse Fourier Transforms of this size
	double* realspace; //< array for holding real-space side of transform data
	complex<double>* kspace; //< array for holding kspace-side of transform data
};

/// Circulant matrices
/** A circulant matrix is a square matrix in which each row is a cyclic permutation by one of the previous row, e.g.
 \f$ \left| \begin{array}{ccc} a & b & c \\ c & a & b \\ b & c & a \end{array} \right| \f$.
 These matrices are convolution operators on vectors, thus they commute and are diagonalized by a Fourier transform.
 The CMatrix class transparently handles converting circulant matrices into and out of the
 Fourier basis, allowing for computationally efficient handling of matrix operations
 (multiplication, inversion, etc.) of circulant matrices.
 Note, the internal data representation is the transpose of the matrix as defined above;
 	the necessary permutation of component order is automatically applied for vector multiplication.
 The FFTs are performed by the <a href="http://www.fftw.org">FFTW library</a>,
 which pre-calculates plans to expedite FFT'ing specific length data arrays. The CMatrix class keeps a cache of
 the FFTW data needed for each size of CMatrix instantiated (which could become inefficient if a wide variety of
 CMatrix sizes are used in the same code, but is suitable for this application where only one size of CMatrix is used
 for a particular interaction symmetry).*/
class CMatrix {
public:
	/// Constructor
	CMatrix(unsigned int ncyc = 0);
	
	/// Make a CMatrix using a data array for the first column
	static CMatrix cmatrixFromColumn(int ncyc, double* coldat);
	
	/// Save matrix to a file (to be read by readFromFile())
	void writeToFile(std::ostream& o) const;
	/// Read matrix from a file written by writeToFile()
	static CMatrix readFromFile(std::istream& s);
	
	/// generate an identity CMatrix
	static CMatrix identity(unsigned int nc);
	/// Fill the first row of this matrix with the ascending sequence \f$ r_0,r_0+1,r_0+2,\cdots \f$
	static CMatrix ramp(unsigned int nc, double r0);
	/// Fill this CMatrix with random numbers in [0,1]
	static CMatrix random(unsigned int nc);
	
	/// zero all entries in this CMatrix
	void zero();
	/// Print this CMatrix to stdout
	void display() const;
	/// Print kspace data for this CMatrix to stdout
	void displayK() const;
	
	/// immutable element access
	double operator[](unsigned int n) const;
	/// mutable element access
	double& operator[](unsigned int n);
	
	/// L2 (Spectral) norm of circulant matrix
	double norm_L2() const;
	/// determinant of circulant matrix
	double det() const;
	/// trace of circulant matrix
	double trace() const;
	
	/// Return a pointer to the CMatrix's Fourier representation
	std::vector< complex<double> >& getKData();
	/// Return a pointer to the CMatrix's Fourier representation (read only)
	const std::vector< complex<double> >& getKData() const;
	/// Return a pointer to the CMatrix's real-space representation
	std::vector<double>& getRealData();
	/// Return a pointer to the CMatrix's real-space representation (read only)
	const std::vector<double>& getRealData() const;
			
	/// Calculate the inverse of this CMatrix
	const CMatrix inverse() const;
	/// Invert this CMatrix inplace
	CMatrix& invert();
	/// Return the transpose of this CMatrix
	const CMatrix transpose() const;
	
	/// unary minus
	const CMatrix operator-() const;
	
	/// Product with a scalar, inplace
	CMatrix& operator*=(double c);
	/// Product with another CMatrix (circulant matrices are commutative!)
	CMatrix& operator*=(const CMatrix& m);
	/// Product with a scalar
	const CMatrix operator*(double c) const;
	/// Product with another CMatrix (circulant matrices are commutative!)
	const CMatrix operator*(const CMatrix& m) const;
	/// Multiply a vector on the right
	const VarVec<double> operator*(const VarVec<double>& v) const;
	
	/// add another CMatrix to this one
	CMatrix& operator+=(const CMatrix& rhs);
	/// sum of two CMatrices
	const CMatrix operator+(const CMatrix& rhs) const;
	/// subtract another CMatrix from this one
	CMatrix& operator-=(const CMatrix& rhs);
	/// difference of two CMatrices
	const CMatrix operator-(const CMatrix& rhs) const;
	
	/// Print the rth row of the matrix to stdout
	void printRow(int r) const;
	
	unsigned int ncycles; //< number of rows (columns)
	
private:
	
	/// Allocate memory for the real-space data of the matrix
	void alloc_data() const;
	/// Allocate memory for the k-space data of the matrix
	void alloc_kdata() const;
	/// create a new FFTW plan for a different size of CMatrix
	void add_inverter() const;
	/// get the appropriate FFTW plans for this size of CMatrix
	cmatrix_fft& get_ffter() const;
	
	/// calculate K-space data from real space
	void calculateKData() const;
	/// calculate real-space data from K-space
	void calculateRealData() const;
	
	mutable std::vector<double> data;				//< real-space data
	mutable std::vector< complex<double>> kdata;	//< K-space data
	mutable bool has_realspace;						//< whether the real-space representation of this matrix has been calculated
	mutable bool has_kspace;						//< whether the k-space representation of this matrix has been calculated
	
	static std::vector<cmatrix_fft> ffters;	//< cache of FFTW plans for FFT'ing various sizes of CMatrix
};

/// string format for CMatrix display
std::ostream& operator<<(std::ostream& o, const CMatrix& m);

#endif
