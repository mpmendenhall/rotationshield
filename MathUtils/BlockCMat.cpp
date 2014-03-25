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
	std::sort(svalues.getData().begin(), svalues.getData().end());
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
		const VarMat<lapack_complex_double>& bPsI = block_SVDs[i]->calc_pseudo_inverse(epsilon*svalues.getData().back());
		for(unsigned int r=0; r<M; r++) {
			for(unsigned int c=0; c<N; c++) {
				(*PsI)(r,c).getKData()[i] = bPsI(r,c);
			}
		}
	}
#endif
	return *PsI;
}

void BlockCMat_SVD::writeToFile(std::ostream& o) const {
	writeString("(BlockCMat_SVD)",o);
	o.write((char*)&M,				sizeof(M));
	o.write((char*)&N,				sizeof(N));
	o.write((char*)&Mc,				sizeof(Mc));
#ifdef WITH_LAPACKE
	assert(block_SVDs.size() == Mc/2+1);
	for(unsigned int i=0; i<block_SVDs.size(); i++)
		block_SVDs[i]->writeToFile(o);
#endif
	svalues.writeToFile(o);
	o.write((char*)&PsI,			sizeof(PsI));
	if(PsI) PsI->writeToFile(o);
	o.write((char*)&PsI_epsilon,	sizeof(PsI_epsilon));
	writeString("(/BlockCMat_SVD)",o);
}

BlockCMat_SVD* BlockCMat_SVD::readFromFile(std::istream& s) {
	checkString("(BlockCMat_SVD)",s);
	BlockCMat_SVD* foo = new BlockCMat_SVD();
	s.read((char*)&foo->M,			sizeof(foo->M));
	s.read((char*)&foo->N,			sizeof(foo->N));
	s.read((char*)&foo->Mc,			sizeof(foo->Mc));
#ifdef WITH_LAPACKE
	for(unsigned int i=0; i<foo->Mc/2+1; i++)
		foo->block_SVDs.push_back( LAPACKE_Matrix_SVD<double,lapack_complex_double>::readFromFile(s) );
#endif
	foo->svalues = VarVec<double>::readFromFile(s);
	s.read((char*)&foo->PsI,		sizeof(foo->PsI));
	if(foo->PsI) {
		foo->PsI = new BlockCMat;
		*foo->PsI = BlockCMat::readFromFile(s);
	}
	s.read((char*)&foo->PsI_epsilon,sizeof(foo->PsI_epsilon));
	checkString("(/BlockCMat_SVD)",s);
	return foo;
}
