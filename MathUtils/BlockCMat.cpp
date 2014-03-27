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

// utility class for sorting enumerated singular values
class Compare_BCM_SVD_singular_values {
public:
	Compare_BCM_SVD_singular_values(const BlockCMat_SVD* b): B(b) {}
	bool operator() (unsigned int i, unsigned int j) { return B->getSV(i) < B->getSV(j); }
	const BlockCMat_SVD* B;
};


BlockCMat_SVD::BlockCMat_SVD(const BlockCMat& BC): M(BC.nRows()), N(BC.nCols()), Mc(BC[0].nRows()), Ms(std::min(M,N)), PsI(NULL), PsI_epsilon(0) {
#ifdef WITH_LAPACKE
	for(unsigned int i=0; i<Mc/2+1; i++) {
		VarMat<lapack_complex_double> dblock(BC.nRows(),BC.nCols());
		for(unsigned int r=0; r<M; r++) {
			for(unsigned int c=0; c<N; c++) {
				dblock(r,c) = BC(r,c).getKData()[i];
			}
		}
		block_SVDs.push_back(new LAPACKE_Matrix_SVD<double,lapack_complex_double>(dblock));
	}
	sort_singular_values();
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

double BlockCMat_SVD::getSV(unsigned int i) const {
#ifdef WITH_LAPACKE
	unsigned int idiag = i/Ms;
	assert(idiag < Mc);
	if(idiag>=Mc/2+1) idiag = Mc - idiag;
	return block_SVDs[idiag]->singular_values()[i%Ms];
#else
	assert(false);
	return 0;
#endif
}

void BlockCMat_SVD::sort_singular_values() {
	// clear any previous sorted lists
	svalues.getData().clear();
	sloc.getData().clear();
	
	// sort internal singular value "ID numbers"
	for(unsigned int i=0; i<Ms*Mc; i++) sloc.push_back(i);
	Compare_BCM_SVD_singular_values CB(this);
	std::sort(sloc.getData().begin(), sloc.getData().end(), CB);
	
	// sorted list of singular values
	for(unsigned int i=0; i<Ms*Mc; i++) svalues.push_back(getSV(sloc[i]));
}

VarVec<double> BlockCMat_SVD::getRightSVec(unsigned int i) const {
	i = sloc[i];
	VarVec<double> v;
#ifdef WITH_LAPACKE
	// check whether this vector belongs to implied complementary set
	unsigned int idiag = i/Ms;
	assert(idiag < Mc);
	bool iset = (idiag>Mc/2+1) || (idiag==Mc/2+1 && !(Mc%2));
	if(iset) idiag = Mc - idiag;
	
	VarVec<lapack_complex_double> sv = block_SVDs[idiag]->getRightSVec(i%Ms);
	assert(sv.size()==M);
	for(unsigned int m=0; m<M; m++) {
		CMatrix C(Mc);
		C.getKData()[idiag] = sv[m] * (iset ? complex<double>(0,1) : 1);
		const std::vector<double>& r = C.getRealData();
		for(unsigned int mc=0; mc<Mc; mc++) v.push_back(r[mc]);
	}
#endif
	return v;
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
	foo->Ms = std::min(foo->M,foo->N);
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
	foo->sort_singular_values();
	return foo;
}
