/// \file "CMatrix.hh" \brief Circulant matrices
#ifndef CMATRIX_HH
/// Make sure this header is only loaded once
#define CMATRIX_HH 1

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include "VarVec.hh"
#include "ComplexT.hh"

typedef ComplexT<double> cdouble;

/// Stores fftw data for FFT'ing 
struct cmatrix_fft
{
	unsigned int ncyc; //< size of CMatrix this is meant for
	fftw_plan forwardplan; //< FFTW data for forward Fourier Transforms of this size
	fftw_plan reverseplan; //< FFTW data for inverse Fourier Transforms of this size
	double* realspace; //< array for holding real-space side of transform data
	cdouble* kspace; //< array for holding kspace-side of transform data
};

/// Circulant matrices
/** A circulant matrix is a square matrix in which each row is a cyclic permutation by one of the previous row, e.g.
 \f$ \left| \begin{array}{ccc} a & b & c \\ b & c & a \\ c & a & b \end{array} \right| \f$.
 These matrices are convolution operators on vectors, thus they commute and are diagonalized by a Fourier transform.
 The CMatrix class transparently handles converting circulant matrices into and out of the
 Fourier basis, allowing for computationally efficient handling of matrix operations
 (multiplication, inversion, etc.) of circulant matrices. The FFTs are performed by the <a href="http://www.fftw.org">FFTW library</a>,
 which pre-calculates plans to expedite FFT'ing specific length data arrays. The CMatrix class keeps a cache of
 the FFTW data needed for each size of CMatrix instantiated (which could become inefficient if a wide variety of
 CMatrix sizes are used in the same code, but is suitable for this application where only one size of CMatrix is used
 for a particular shield simulation).*/
template<typename T>
class CMatrix {
public:
	/// Constructor
	CMatrix(unsigned int ncyc = 0);
	/// COnstructor from data
	template <class InputIterator>
	CMatrix(InputIterator first, InputIterator last);
	/// Destructor
	~CMatrix() {}
	
	/// Save matrix to a file (to be read by readFromFile())
	void writeToFile(std::ostream& o) const;
	/// Read matrix from a file written by writeToFile()
	static CMatrix<T> readFromFile(std::istream& s);
	
	/// generate an identity CMatrix
	static CMatrix<T> identity(unsigned int nc);
	/// Fill the first row of this matrix with the ascending sequence \f$ r_0,r_0+1,r_0+2,\cdots \f$
	static CMatrix<T> ramp(unsigned int nc, T r0);
	/// Fill this CMatrix with random numbers in [0,1]
	static CMatrix<T> random(unsigned int nc);
	
	/// zero all entries in this CMatrix
	void zero();
	/// Print this CMatrix to stdout
	void display() const;
	/// Print kspace data for this CMatrix to stdout
	void displayK() const;
	
	/// immutable element access
	T operator[](unsigned int n) const;
	/// mutable element access
	T& operator[](unsigned int n);
	
	/// L2 (Spectral) norm of circulant matrix
	T norm_L2() const;
	
	/// Return a pointer to the CMatrix's Fourier representation
	std::vector<cdouble>& getKData();
	/// Return a pointer to the CMatrix's Fourier representation (read only)
	const std::vector<cdouble>& getKData() const;	
	/// Return a pointer to the CMatrix's real-space representation
	std::vector<T>& getRealData();
	/// Return a pointer to the CMatrix's real-space representation (read only)
	const std::vector<T>& getRealData() const;
		
	/// Allocate memory for the real-space data of the matrix
	void alloc_data() const;
	/// Allocate memory for the k-space data of the matrix
	void alloc_kdata() const;
	
	/// Calculate the inverse of this CMatrix
	const CMatrix<T> inverse() const;
	/// Invert this CMatrix inplace
	CMatrix<T>& invert();
	/// Return the transpose of this CMatrix
	const CMatrix<T> transpose() const;
	
	/// unary minus
	const CMatrix<T> operator-() const;
	
	/// Product with a scalar, inplace
	CMatrix<T>& operator*=(T c);
	/// Product with another CMatrix (circulant matrices are commutative!)
	CMatrix<T>& operator*=(const CMatrix<T>& m);
	/// Product with a scalar
	const CMatrix<T> operator*(T c) const;
	/// Product with another CMatrix (circulant matrices are commutative!)
	const CMatrix<T> operator*(const CMatrix<T>& m) const;
	/// Multiply a vector on the right
	const VarVec<T> operator*(const VarVec<T>& v) const;
	
	/// add another CMatrix to this one
	CMatrix<T>& operator+=(const CMatrix& rhs);
	/// sum of two CMatrices
	const CMatrix<T> operator+(const CMatrix& rhs) const;
	/// subtract another CMatrix from this one
	CMatrix<T>& operator-=(const CMatrix& rhs);
	/// difference of two CMatrices
	const CMatrix<T> operator-(const CMatrix& rhs) const;
	
	
	/// Print the rth row of the matrix to stdout
	void printRow(int r) const;
	/// Make a CMatrix using a data array for the first column
	static CMatrix<T> cmatrixFromColumn(int ncyc, T* coldat);
	
	unsigned int ncycles; //< number of rows (columns)
	
	void calculateKData() const;
	void calculateRealData() const;
	
private:
	
	mutable std::vector<T> data;
	mutable std::vector<cdouble> kdata;
	mutable bool has_realspace; //< whether the real-space representation of this matrix has been calculated
	mutable bool has_kspace; //< whether the k-space representation of this matrix has been calculated
	static std::vector<cmatrix_fft> ffters; //< cache of FFTW plans for FFT'ing various sizes of CMatrix
	void add_inverter() const; //< create a new FFTW plan for a different size of CMatrix
	cmatrix_fft& get_ffter() const; //< get the appropriate FFTW plans for this size of CMatrix
};

template<typename T>
std::vector<cmatrix_fft> CMatrix<T>::ffters = std::vector<cmatrix_fft>();

// Constructor and Destructor ------______------______-------_______------______------

template<typename T>
CMatrix<T>::CMatrix(unsigned int ncyc): ncycles(ncyc), data(), kdata(), has_realspace(false), has_kspace(false) {
}

template<typename T>
template <class InputIterator>
CMatrix<T>::CMatrix(InputIterator first, InputIterator last): ncycles(last-first), data(first,last), kdata(), has_realspace(true), has_kspace(false) {
}

template<typename T>
void CMatrix<T>::writeToFile(std::ostream& o) const {
	o.write((char*)&ncycles,sizeof(int));
	o.write((char*)&has_realspace,sizeof(bool));
	o.write((char*)&has_kspace,sizeof(bool));
	if(has_realspace)
		o.write((char*)&data.front(),sizeof(T)*ncycles);
	if(has_kspace)
		o.write((char*)&kdata.front(),sizeof(cdouble)*(ncycles/2+1));
}

template<typename T>
CMatrix<T> CMatrix<T>::readFromFile(std::istream& s) {
	int ncyc;
	bool hrd, hkd;
	s.read((char*)&ncyc,sizeof(int));
	s.read((char*)&hrd,sizeof(bool));
	s.read((char*)&hkd,sizeof(bool));
	CMatrix<T> foo = CMatrix<T>(ncyc);
	if(hrd) {
		foo.alloc_data();
		s.read((char*)&foo.data.front(),sizeof(T)*ncyc);
	}
	if(hkd) {
		foo.alloc_kdata();
		s.read((char*)&foo.kdata.front(),sizeof(cdouble)*(ncyc/2+1));
	}
	return foo;
}

// Special Matrices ------______------______-------_______------______------

template<typename T>
CMatrix<T> CMatrix<T>::identity(unsigned int nc) {
	CMatrix<T> m = CMatrix<T>(nc);
	m[0] = 1.0;
	return m;
}

template<typename T>
void CMatrix<T>::zero() {
	has_kspace = true;
	has_realspace = true;
	data = std::vector<T>(ncycles);
	kdata = std::vector<cdouble>(ncycles/2+1);
}

template<typename T>
CMatrix<T> CMatrix<T>::ramp(unsigned int nc, T r0) {
	CMatrix<T> m = CMatrix<T>(nc);
	for(unsigned int i=0; i<nc; i++)
		m[i] = i+r0;
	return m;
}

template<typename T>
CMatrix<T> CMatrix<T>::random(unsigned int nc) {
	CMatrix<T> m = CMatrix<T>(nc);
	for(unsigned int i=0; i<nc; i++)
		m[i] = T(rand())/T(RAND_MAX);
	return m;
}

// Data Providers ------______------______-------_______------______------

template<typename T>
void CMatrix<T>::alloc_data() const {
	data = std::vector<T>(ncycles);
	has_realspace = true;
}

template<typename T>
void CMatrix<T>::alloc_kdata() const {
	kdata = std::vector<cdouble>(ncycles);
	has_kspace = true;
}


template<typename T>
void CMatrix<T>::calculateRealData() const {
	cmatrix_fft& f = get_ffter();
	for(unsigned int i=0; i<ncycles/2+1; i++)
		f.kspace[i] = kdata[i];
	fftw_execute(f.reverseplan);
	alloc_data();
	for(unsigned int n=0; n<ncycles; n++) data[n] = T(f.realspace[n])/T(ncycles);	
}

template<typename T>
T CMatrix<T>::operator[](unsigned int i) const {
	assert(i<ncycles);
	if(has_realspace)
		return data[i];
	if(!has_kspace) {
		alloc_data();
		return data[i];
	}
	calculateRealData();
	return data[i];
}

template<typename T>
T& CMatrix<T>::operator[](unsigned int i) {
	assert(i<ncycles);
	if(has_realspace)
		return data[i];
	if(!has_kspace) {
		alloc_data();
		return data[i];
	}
	calculateRealData();
	has_kspace = false; // messing with real space data invalidates kspace data
	return data[i];	
}


template<typename T>
void CMatrix<T>::calculateKData() const {
	cmatrix_fft& f = get_ffter();
	for(unsigned int n=0; n<ncycles; n++) f.realspace[n] = (double)data[n];
	fftw_execute(f.forwardplan);
	alloc_kdata();
	kdata.assign(f.kspace,f.kspace+ncycles/2+1);
	has_kspace = true;
}

template<typename T>
std::vector<cdouble>& CMatrix<T>::getKData() {
	if(has_kspace)
		return kdata;
	if(!has_realspace) {
		alloc_kdata();
		return kdata;
	}
	calculateKData();
	has_realspace = false;
	return kdata;
}

template<typename T>
const std::vector<cdouble>& CMatrix<T>::getKData() const {
	if(has_kspace)
		return kdata;
	if(!has_realspace) {
		alloc_kdata();
		return kdata;
	}
	calculateKData();
	return kdata;
}

template<typename T>
std::vector<T>& CMatrix<T>::getRealData() {
	if(has_realspace)
		return data;
	if(!has_kspace) {
		alloc_data();
		return data;
	}
	calculateRealData();
	has_kspace = false;
	return data;
}

template<typename T>
const std::vector<T>& CMatrix<T>::getRealData() const {
	if(has_realspace)
		return data;
	if(!has_kspace) {
		alloc_data();
		return data;
	}
	calculateRealData();
	return data;
}

// Matrix properties ------______------______-------_______------______------

template<typename T>
T CMatrix<T>::norm_L2() const {
	const std::vector<cdouble>& v = getKData();
	std::vector<double> vn;
	for(std::vector<cdouble>::const_iterator it = v.begin(); it < v.end(); it++)
		vn.push_back(it->mag());
	return *std::max_element(vn.begin(),vn.end());
}


// "GUI" ------______------______-------_______------______------

template<typename T>
void CMatrix<T>::printRow(int r) const {
	printf("| ");
	for(unsigned int i=0; i<ncycles; i++) printf("%.3g ",(double)(*this)[(i+(ncycles-r))%ncycles]);
	printf("|");
}

template<typename T>
void CMatrix<T>::display() const {
	std::cout << "CMatrix " << ncycles << " " << has_realspace << " " << has_kspace << std::endl;
	for(int r=0; r<ncycles; r++) {
		printRow(r);
		printf("\n");
	}
	std::cout << "CMatrix " << ncycles << " " << has_realspace << " " << has_kspace << std::endl;
}

template<typename T>
void CMatrix<T>::displayK() const {
	std::cout << "{ ";
	for(int i=0; i<ncycles/2+1; i++) std::cout << getKData()[i] << " ";
	std::cout << "}" << std::endl;
}

/// string format for CMatrix display
template<typename T>
std::ostream& operator<<(std::ostream& o, const CMatrix<T>& m) {
	for(unsigned int r=0; r<m.ncycles; r++) {
		o << "| ";
		for(int c=0; c<m.ncycles; c++)
			o << m[(c+(m.ncycles-r))%m.ncycles] << " ";
		o << "|\n";
	}
	return o;
}

// Matrix Ops ------______------______-------_______------______------

template<typename T>
CMatrix<T>& CMatrix<T>::operator+=(const CMatrix<T>& m)
{
	assert(ncycles == m.ncycles);
	bool kdata_ok, realdata_ok;
	kdata_ok = realdata_ok = false;
	
	if(!has_realspace && !has_kspace) {
		*this = m;
		return *this;
	}
	
	if(has_kspace && (m.has_kspace || !has_realspace)) {
		const std::vector<cdouble>& kd = m.getKData();
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] += kd[i];
		kdata_ok = true;
		has_kspace = true;
	}

	if(has_realspace && (m.has_realspace || !has_kspace)) {
		for(unsigned int i=0; i<ncycles; i++)
			data[i] += m[i];
		realdata_ok = true;
	}


	
	has_realspace = realdata_ok;
	has_kspace = kdata_ok;
	
	return *this;
}

template<typename T>
const CMatrix<T> CMatrix<T>::operator+(const CMatrix<T>& m) const {
	CMatrix<T> r = *this;
	r += m;
	return r;
}

template<typename T>
CMatrix<T>& CMatrix<T>::operator-=(const CMatrix<T>& m)
{
	assert(ncycles == m.ncycles);
	bool kdata_ok, realdata_ok;
	kdata_ok = realdata_ok = false;
	
	if(!has_realspace && !has_kspace) {
		*this = -m;
		return *this;
	}
	
	if(has_kspace && (m.has_kspace || !has_realspace)) {
		const std::vector<cdouble>& kd = m.getKData();
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] -= kd[i];
		kdata_ok = true;
		has_kspace = true;
	}

	if(has_realspace && (m.has_realspace || !has_kspace)) {
		for(unsigned int i=0; i<ncycles; i++)
			data[i] -= m[i];
		realdata_ok = true;
	}
	
	
	has_realspace = realdata_ok;
	has_kspace = kdata_ok;
	
	return *this;
}

template<typename T>
const CMatrix<T> CMatrix<T>::operator-(const CMatrix<T>& m) const {
	CMatrix<T> r = *this;
	r -= m;
	return r;
}

template<typename T>
CMatrix<T>& CMatrix<T>::operator*=(T c) {
	
	if(has_realspace)
		for(unsigned int i=0; i<ncycles; i++) data[i] *= c;
	
	if(has_kspace)
		for(unsigned int i=0; i<ncycles/2+1; i++)
			kdata[i] *= (double)c;
	
	return *this;
}

template<typename T>
const CMatrix<T> CMatrix<T>::operator*(T c) const {
	CMatrix<T> r = *this;
	r *= c;
	return r;
}

template<typename T>
const CMatrix<T> CMatrix<T>::operator-() const {
	CMatrix<T> r = *this;
	r *= -1.0;
	return r;
}

template<typename T>
CMatrix<T>& CMatrix<T>::operator*=(const CMatrix<T>& m) {
	getKData();
	const std::vector<cdouble>& mkd = m.getKData();
	has_realspace = false;
	for(unsigned int i=0; i<ncycles/2+1; i++)
		kdata[i] *= mkd[i];
	return *this;
}

template<typename T>
const CMatrix<T> CMatrix<T>::operator*(const CMatrix<T>& m) const {
	CMatrix<T> r = *this;
	r *= m;
	return r;
}


template<typename T>
const VarVec<T> CMatrix<T>::operator*(const VarVec<T>& v) const
{
	
	const std::vector<cdouble>& kd = getKData();
	
	cmatrix_fft& f = get_ffter();
	f.realspace[0] = (double)v[0];
	for(unsigned int i=1; i<ncycles; i++) f.realspace[i] = (double)v[ncycles - i];
	fftw_execute(f.forwardplan);
	
	for(unsigned int i=0; i<ncycles/2+1; i++)
		f.kspace[i] *= kd[i];
	
	fftw_execute(f.reverseplan);
	
	VarVec<T> out = VarVec<T>(ncycles);
	out[0] = (T)f.realspace[0];
	for(unsigned int i=1; i<ncycles; i++) out[i] = (T)f.realspace[ncycles-i];
	out /= T(ncycles);
	return out;
}

template<typename T>
CMatrix<T>& CMatrix<T>::invert() {
	std::vector<cdouble>& kd = getKData();
	for(unsigned int i=0; i<ncycles/2+1; i++)
		kd[i] = kd[i].inverse();
	return *this;
}

template<typename T>
const CMatrix<T> CMatrix<T>::inverse() const
{
	CMatrix<T> I = *this;
	I.invert();
	return I;
}

template<typename T>
CMatrix<T> CMatrix<T>::cmatrixFromColumn(int ncyc, T* coldat)
{
	CMatrix<T> m = CMatrix<T>(ncyc);
	m[0] = coldat[0];
	for(int n=1; n<ncyc; n++) m[n] = coldat[ncyc-n];
	return m;
}

template<typename T>
const CMatrix<T> CMatrix<T>::transpose() const
{
	CMatrix<T> m = CMatrix<T>(ncycles);
	std::vector<T>& md = m.getRealData(); 
	md[0] = (*this)[0];
	for(int n=1; n<ncycles; n++) md[n] = (*this)[ncycles-n];
	return m;
}

// Inverter handler ------______------______-------_______------______------

template<typename T>
void CMatrix<T>::add_inverter() const
{
	cmatrix_fft NI;
	ffters.push_back(NI);
	cmatrix_fft& newInverter = ffters.back();
	
	newInverter.ncyc = ncycles;
	newInverter.realspace = new double[ncycles];
	newInverter.kspace = new cdouble[ncycles/2+1];
	
	FILE* fin = fopen("fftw_wisdom","r");
	if(fin) {
		fftw_import_wisdom_from_file(fin);
		fclose(fin);
	}
	
	newInverter.forwardplan = fftw_plan_dft_r2c_1d(ncycles,
												   newInverter.realspace,
												   (fftw_complex*)newInverter.kspace,
												   FFTW_EXHAUSTIVE);
	newInverter.reverseplan = fftw_plan_dft_c2r_1d(ncycles,
												   (fftw_complex*)newInverter.kspace,
												   newInverter.realspace,
												   FFTW_EXHAUSTIVE);
	
	FILE* fout = fopen("fftw_wisdom","w");
	if(fout) {
		fftw_export_wisdom_to_file(fout);
		fclose(fout);
	}
}

template<typename T>
cmatrix_fft& CMatrix<T>::get_ffter() const
{
	std::vector<cmatrix_fft>::iterator it;
	for(it = ffters.begin(); it != ffters.end(); it++)
		if(it->ncyc == ncycles) return *it;
	add_inverter();
	return ffters[ffters.size()-1];
}

#endif
