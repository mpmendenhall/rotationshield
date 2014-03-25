#include "BlockCMat.hh"

BlockCMat makeBlockCMatIdentity(unsigned int n, unsigned int mc) {
	BlockCMat foo(n,n,mc);
	for(unsigned int i=0; i<n; i++)
		foo(i,i) = CMatrix::identity(mc);
	return foo;
}

BlockCMat makeBlockCMatRandom(unsigned int n, unsigned int mc) {
	BlockCMat foo(n,n,mc);
	for(unsigned int r=0; r<n; r++)
		for(unsigned int c=0; c<n; c++)
			foo(r,c) = CMatrix::random(mc);
	return foo;
}

//-----------------------------------

BlockCMat_SVD::BlockCMat_SVD(const BlockCMat& BC): M(BC.nRows()), N(BC.nCols()), Mc(BC[0].nRows()), PsI(NULL), PsI_epsilon(0) {
#ifdef WITH_LAPACKE
	for(unsigned int i=0; i<Mc/2+1; i++) {
		VarMat<lapack_complex_double> dblock(BC.nRows(),BC.nCols());
		for(unsigned int r=0; r<M; r++) {
			for(unsigned int c=0; c<N; c++) {
				dblock(r,c) = BC(r,c).getKData()[i];
			}
		}
		block_SVDs.push_back(new LAPACKE_Matrix_SVD<double,lapack_complex_double>(dblock));
		const VarMat<double> S = block_SVDs.back()->singular_values();
		for(unsigned int j=0; j<S.size(); j++) svalues.push_back(S[j]);
	}
	std::sort(svalues.begin(),svalues.end());
#else
	PsI = new BlockCMat(BC);
	PsI->invert();
#endif
}

BlockCMat_SVD::~BlockCMat_SVD() {
	if(PsI) delete PsI;
#ifdef WITH_LAPACKE
	for(unsigned int i=0; i<block_SVDs.size(); i++) delete(block_SVDs[i]);
#endif
}
	
const BlockCMat& BlockCMat_SVD::calc_pseudo_inverse(double epsilon) {
#ifdef WITH_LAPACKE
	if(PsI && PsI_epsilon==epsilon) return *PsI;
	if(PsI) delete PsI;
	PsI_epsilon = epsilon;
	
	PsI = new BlockCMat(M, N, Mc);
	for(unsigned int i=0; i<Mc/2+1; i++) {
		const VarMat<lapack_complex_double>& bPsI = block_SVDs[i]->calc_pseudo_inverse(epsilon*svalues.back());
		for(unsigned int r=0; r<M; r++) {
			for(unsigned int c=0; c<N; c++) {
				(*PsI)(r,c).getKData()[i] = bPsI(r,c);
			}
		}
	}
#endif
	return *PsI;
}
