/// \file "BlockCMat.hh" \brief Block Circulant matrices
#ifndef BLOCKCMAT_HH
/// Make sure this header is only loaded once
#define BLOCKCMAT_HH 1

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cassert>
#include "CMatrix.hh"
#include "VarVec.hh"

/// Block Circulant matrix: Deprecated! Use VarMat< CMatrix<T> >
/** A Block Circulant matrix is a matrix which can be subdivided into equal-size blocks which are circulant matrices
 (stored as a CMatrix). Since the circulant matrix components are commutative, a \f$ (n\cdot c)\times (m\cdot c)\f$
 block circulant matrix composed of \f$ n \times m \f$ \f$ c \times c \f$ circulant blocks can be handled much like
 a \f$ n \times m \f$ matrix of scalars. Computations with block circulant matrices are faster than those with
 general matrices of the same dimensions beause of the faster handling of computations with the circulant submatrices. */
template<typename T>
class BlockCMat
	{
	public:
		/// Constructor \param nrw number of rows of circulants \param ncl number of columns of circulants \param ncy dimension of the circulant blocks
		BlockCMat(int nrw = 0, int ncl = 0, int ncy = 0);
		/// Destructor
		~BlockCMat() {}
		
		/// Write data to a file (to be read by readFromFile())
		void writeToFile(std::ostream& o) const;
		/// Read a BlockCMat from a file written by writeToFile()
		static BlockCMat readFromFile(std::istream& s);
		/// generate an Identity block-circulant matrix
		static BlockCMat identity(unsigned int nrw, unsigned int ncy);
		/// Fill this matrix with random values in [0,1]
		static BlockCMat random(unsigned int nrw, unsigned int ncl, unsigned int ncy);
		
		
		/// Print this matrix to stdout
		void display() const;
		/// Print the rth row to stdout
		void printRow(unsigned int r) const;
		/// Print info about this matrix to stdout
		void dispInfo() const { printf("block matrix: %i x %i %i x %i circulants\n",nrows,ncols,ncycles,ncycles); }
		
		///unary minus
		const BlockCMat<T> operator-() const;
		/// matrix product
		const BlockCMat<T> operator*(const BlockCMat<T>& rhs) const;
		/// matrix sum
		const BlockCMat<T> operator+(const BlockCMat<T>& rhs) const;
		/// matrix difference
		const BlockCMat<T> operator-(const BlockCMat<T>& rhs) const;
		/// matrix right-multiplied by a vector
		const VarVec<T> operator*(const VarVec<T>& v) const;
		
		/// Get the specified circulant block
		CMatrix<T>& getBlock(unsigned int r, unsigned int c) { assert(r<nrows && c<ncols); return data[r*ncols + c]; }
		/// Get the specified circulant block, read only
		const CMatrix<T>& getBlock(unsigned int r, unsigned int c) const { assert(r<nrows && c<ncols); return data[r*ncols+c]; }
		/// Get a submatrix of this matrix
		const BlockCMat<T> getSubmatrix(unsigned int r0, unsigned int r1, unsigned int c0, unsigned int c1) const;
		
		/// Check how well two BlockCMats are inverses of each other
		static T checkInversion(BlockCMat<T>& a, BlockCMat<T>& b);
		
		/// Invert this matrix, inplace
		void invert();
		
	private:
		
		unsigned int nrows; //< number of rows of circulants
		unsigned int ncols; //< number of columns of circulants
		unsigned int ncycles; //< dimensions of the circulant blocks
		std::vector< CMatrix<T> > data; //< the circulant blocks
		
		/// Part of the inversion process
		void invert_firstcell_inplace();
		/// Part of the inversion process
		void clear_firstcolumn_inplace();
	};

//------______------______-------_______------______------

template<typename T>
BlockCMat<T>::BlockCMat(int nrw,int ncl, int ncy): nrows(nrw), ncols(ncl), ncycles(ncy)  {
	data = std::vector< CMatrix<T> >(nrows*ncols,CMatrix<T>(ncycles));
}

template<typename T>
const BlockCMat<T> BlockCMat<T>::getSubmatrix(unsigned int r0,unsigned int r1, unsigned int c0, unsigned int c1) const {
	assert(r0 <= r1 && r1 < nrows && c0 <= c1 && c1 < ncols);
	BlockCMat<T> s = BlockCMat<T>(1+r1-r0, 1+c1-c0, ncycles);
	for(unsigned int r=0; r<s.nrows; r++)
		for(unsigned int c=0; c<s.ncols; c++)
			s.getBlock(r,c) = getBlock(r0+r,c0+c);
	return s;
}

template<typename T>
void BlockCMat<T>::writeToFile(std::ostream& o) const {
	o.write((char*)&nrows,sizeof(int));
	o.write((char*)&ncols,sizeof(int));
	o.write((char*)&ncycles,sizeof(int));	
	for(unsigned int r=0; r<nrows; r++)
		for(unsigned int c=0; c<ncols; c++)
			getBlock(r,c).writeToFile(o);
}

template<typename T>
BlockCMat<T> BlockCMat<T>::readFromFile(std::istream& s)
{
	int nr, nc, ncyc;
	s.read((char*)&nr,sizeof(int));
	s.read((char*)&nc,sizeof(int));
	s.read((char*)&ncyc,sizeof(int));
	
	BlockCMat<T> foo = BlockCMat<T>(nr,nc,ncyc);
	for(int r=0; r<nr; r++)
		for(int c=0; c<nc; c++)
			foo.getBlock(r,c) = CMatrix<T>::readFromFile(s);
	return foo;
}


template<typename T>
BlockCMat<T> BlockCMat<T>::identity(unsigned int nrw, unsigned int ncy) {
	BlockCMat<T> I = BlockCMat<T>(nrw,nrw,ncy);
	for(int r=0; r<nrw; r++) I.getBlock(r,r) = CMatrix<T>::identity(ncy);
	return I;
}

template<typename T>
BlockCMat<T> BlockCMat<T>::random(unsigned int nrw, unsigned int ncl, unsigned int ncy) {
	BlockCMat<T> R = BlockCMat<T>(nrw,ncl,ncy);
	for(int r=0; r<nrw; r++)
		for(int c=0; c<ncl; c++)
			R.getBlock(r,c) = CMatrix<T>::random(ncy);
	return R;
}


template<typename T>
void BlockCMat<T>::printRow(unsigned int r) const
{
	for(unsigned int r0 = 0; r0<ncycles; r0++)
	{
		printf("|");
		for(unsigned int c=0; c<ncols; c++) getBlock(r,c).printRow(r0);
		printf("|\n");
	}
}

template<typename T>
void BlockCMat<T>::display() const
{
	for(unsigned int r=0; r<nrows; r++) printRow(r);
	printf("\n");
}

template<typename T>
const BlockCMat<T> BlockCMat<T>::operator-() const {
	BlockCMat<T> r = *this;
	for(unsigned int r0=0; r0<nrows; r0++)
		for(unsigned int c0=0; c0<ncols; c0++)
			r.getBlock(r0,c0) *= -1.0;
	return r;
}

template<typename T>
const BlockCMat<T> BlockCMat<T>::operator*(const BlockCMat<T>& m) const {
	assert(ncols == m.nrows);
	BlockCMat<T> r = BlockCMat<T>(nrows,m.ncols,ncycles);
	for(unsigned int r0=0; r0<nrows; r0++) {
		for(unsigned int c0=0; c0<m.ncols; c0++) {
			for(unsigned int i=0; i<ncols; i++) {
				CMatrix<T> foo = getBlock(r0,i) * m.getBlock(i,c0);
				r.getBlock(r0,c0) += foo;
			}
		}
	}
	return r;
}

template<typename T>
const BlockCMat<T> BlockCMat<T>::operator+(const BlockCMat<T>& m) const {
	assert(ncols == m.ncols && nrows == m.nrows);
	BlockCMat<T> foo = *this;
	for(unsigned int r=0; r<nrows; r++)
		for(unsigned int c=0; c<ncols; c++)
			foo.getBlock(r,c) += m.getBlock(r,c);
	return foo;
}

template<typename T>
const BlockCMat<T> BlockCMat<T>::operator-(const BlockCMat<T>& m) const {
	assert(ncols == m.ncols && nrows == m.nrows);
	BlockCMat<T> foo = *this;
	for(unsigned int r=0; r<nrows; r++)
		for(unsigned int c=0; c<ncols; c++)
			foo.getBlock(r,c) -= m.getBlock(r,c);
	return foo;
}

template<typename T>
const VarVec<T> BlockCMat<T>::operator*(const VarVec<T>& v) const
{
	assert(ncols*ncycles == v.size());
	VarVec<T> vnew = VarVec<T>(ncycles*nrows);
	for(unsigned int c=0; c<ncols; c++) {
		VarVec<T> v2 = VarVec<T>(ncycles);
		for(unsigned int i=0; i<ncycles; i++) v2[i] = v[ncycles*c+i];
		for(unsigned int r=0; r<nrows; r++) {
			VarVec<T> v1 = getBlock(r,c) * v2;
			for(unsigned int i=0; i<ncycles; i++) vnew[i+r*ncycles] += v1[i];
		}
	}
	return vnew;
}

template<typename T>
void BlockCMat<T>::invert() 
{
	printf(">"); fflush(stdout);
	
	if(nrows == 1) {
		getBlock(0,0).invert();
		printf("*<"); fflush(stdout);
		return;
	}
	
	
	clear_firstcolumn_inplace();
	
	//invert the submatrix
	BlockCMat<T> submat = getSubmatrix(1,nrows-1,1,ncols-1);
	submat.invert();
	
	BlockCMat<T> firstcol = getSubmatrix(1,nrows-1,0,0);
	BlockCMat<T> colmul = submat * firstcol;
	
	//put back into this matrix
	for(unsigned int c=1; c<ncols; c++)
		for(unsigned int r=1; r<nrows; r++)
			getBlock(r,c) = submat.getBlock(r-1,c-1);
	for(unsigned int r=1; r<nrows; r++)
		getBlock(r,0) = colmul.getBlock(r-1,0);
	
	//finish off by cleaning first row
	BlockCMat<T> toprow = -getSubmatrix(0,0,1,ncols-1);
	
	for(unsigned int c=0; c<ncols; c++) {
		if(c>0)
			getBlock(0,c).zero();
		for(unsigned int r=1; r<nrows; r++)
			getBlock(0,c) += getBlock(r,c) * toprow.getBlock(0,r-1);
	}
	
	printf("<"); fflush(stdout);
}

template<typename T>
void BlockCMat<T>::clear_firstcolumn_inplace() 
{
	invert_firstcell_inplace();
	for(unsigned int r=1; r<nrows; r++)
	{
		CMatrix<T> m0  = -getBlock(r,0);
		for(unsigned int c=1; c<ncols; c++)
		{
			getBlock(r,c) += getBlock(0,c) * m0;
		}
		getBlock(r,0) = m0 * getBlock(0,0);
	}
}

template<typename T>
void BlockCMat<T>::invert_firstcell_inplace() 
{
	getBlock(0,0).invert();
	for(unsigned int i=1; i<ncols; i++)
		getBlock(0,i) *= getBlock(0,0);
}



template<typename T>
T BlockCMat<T>::checkInversion(BlockCMat<T>& a, BlockCMat<T>& b)
{
	BlockCMat<T> ab = a * b;
	T mx = 0.0;
	T tmp;
	for(unsigned int r=0; r<ab.nrows; r++) {
		for(unsigned int c=0; c<ab.ncols; c++) {
			CMatrix<T> i = ab.getBlock(r,c);
			for(unsigned int j=0; j<i.ncycles; j++) {
				if(r == c && j==0) tmp = fabs(i[j]-1.0);
				else tmp = fabs(i[j]);
				if(tmp > mx) mx = tmp;
			}
		}
	}
	return mx;
}

#endif
