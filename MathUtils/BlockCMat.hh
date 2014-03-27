/// \file "BlockCMat.hh" \brief Block Circulant matrices, combining Circulant FFTW with internal LAPACKE matrix calculations
#ifndef BLOCKCMAT_HH
/// Make sure this header is only loaded once
#define BLOCKCMAT_HH 1

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cassert>
#include "CMatrix.hh"
#include "VarMat.hh"
#include "BinaryOutputObject.hh"
#ifdef WITH_LAPACKE
#include "LAPACKE_Matrix.hh"
#endif

typedef VarMat<CMatrix> BlockCMat;

BlockCMat makeBlockCMatIdentity(unsigned int n, unsigned int mc);
BlockCMat makeBlockCMatRandom(unsigned int n, unsigned int mc);

/// singular-value decomposition of block circulant matrix
class BlockCMat_SVD: public BinaryOutputObject {
public:
	/// constructor
	BlockCMat_SVD(const BlockCMat& BC);
	/// destructor
	~BlockCMat_SVD();
	
	/// generate inverse at given singular value threshold
	const BlockCMat& calc_pseudo_inverse(double epsilon = 0);
	
	/// return sorted list of singular values, for threshold determination
	const std::vector<double>& singular_values() const { return svalues.getData(); }
	
	/// Dump binary data to file
	void writeToFile(std::ostream& o) const;
	/// Read binary data from file
	static BlockCMat_SVD* readFromFile(std::istream& s);
	
	/// sub-block singular values
	double getSV(unsigned int i) const;
	/// get enumerated right singular vector
	VarVec<double> getRightSVec(unsigned int i) const;
	
protected:
	/// empty constructor without calculation
	BlockCMat_SVD() {}

	/// prepare sorted singular values
	void sort_singular_values();
	
	unsigned int M, N, Mc, Ms;
#ifdef WITH_LAPACKE
	std::vector< LAPACKE_Matrix_SVD<double,lapack_complex_double>* > block_SVDs;
#endif
	VarVec<double> svalues;		//< sorted singular values
	VarVec<unsigned int> sloc;	//< location of sorted singular-values in sub-matrix
	BlockCMat* PsI;				//< pseudo-inverse
	double PsI_epsilon;			//< threshold for singular vectors
};

#endif
