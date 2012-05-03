/// \file "Permutation.hh" \brief Permutation matrices
#ifndef PERMUTATION_HH
/// Make sure this header is loaded only once
#define PERMUTATION_HH 1

#include <stdlib.h>
#include <vector>
#include <cassert>

/// Permutation matrix (each row and column has one entry which is 1, all others 0)
class Permutation {
public:
	/// Constructor
	Permutation(unsigned int n);
	/// Destructor
	~Permutation() {}
	
	/// get size
	unsigned int size() const { return data.size(); }
	/// immutable element access
	unsigned int operator[](unsigned int n) const { assert(n<size()); return data[n]; }
	/// mutable element access
	unsigned int& operator[](unsigned int n) { assert(n<size()); return data[n]; }
	
	/// Set to permutation matrix needed by ShieldPhysics
	Permutation& nshuffle(int n);
	
	/// Permute two elements
	void swap(unsigned int a, unsigned int b);
	
	/// Return the inverse (=transpose) of the matrix
	Permutation inverse() const;
	/// Return the transpose (=inverse) of the matrix
	Permutation transposed() const { return inverse(); }
	/// Invert (=transpose) the matrix inplace
	Permutation& invert();
	/// Transpose (=invert) the matrix inplace
	Permutation& transpose() {return invert(); }
	
	/// Multiply another Permutation to the right
	const Permutation operator*(const Permutation& p) const;
	
private:
	std::vector<unsigned int> data; //< the permutation data
};

#endif
