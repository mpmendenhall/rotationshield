#include "SymmetricSolver.hh"
#include "ProgressBar.hh"
#include <cassert>
#include <fstream>

void SymmetricSolver::circulantMul(const BlockCMat& M, mvec& v, unsigned int nPhi) {
	assert(!(v.size()%nPhi));
	assert(M.nCols()*nPhi == v.size());
	
	// stuff vector into vector-of-vectors for circulant blocks
	VarVec< mvec > vv;
	for(unsigned int i=0; i<v.size()/nPhi; i++)
		vv.push_back(mvec(&v[i*nPhi], &v[i*nPhi]+nPhi));
	vv = M.lMultiply<mvec,mvec>(vv);
	
	// pull data back out
	v.getData().resize(M.nRows()*nPhi);
	for(unsigned int i=0; i<M.nRows(); i++)
		for(unsigned int j=0; j<nPhi; j++)
			v[i*nPhi+j] = vv[i][j];
}

BlockCMat SymmetricSolver::makeIdentity(unsigned int N, unsigned int nPhi) {
	return BlockCMat::identity(N, CMatrix::identity(nPhi), CMatrix(nPhi));
}

double SymmetricSolver::checkInversion(const BlockCMat& M, const BlockCMat& MI, unsigned int nPhi) {
	return (M * MI - makeIdentity(M.nRows(),nPhi)).getData().max_norm_L2();
}

void SymmetricSolver::solve(ReactiveSet& R) {
	buildInteractionMatrix(R);
	printf("Solving for Green's Function... ");
	BlockCMat M = the_GF; // save a copy for inversion check
	the_GF.invert();
	printf(" %g inversion error.",checkInversion(M,the_GF,R.nPhi));
	printf(" Done.\n");
}

void SymmetricSolver::buildInteractionMatrix(ReactiveSet& R) {
	if(verbose) printf("Building interaction matrix for %i = %i x %i DF...\n", R.nDF(), R.nPhi, R.nDF()/R.nPhi);
	ProgressBar pb = ProgressBar(R.nDF(), R.nPhi, verbose);
	
	the_ixn = BlockCMat(R.nDF()/R.nPhi, R.nDF()/R.nPhi, CMatrix(R.nPhi));
	
	R.startInteractionScan();
	for(unsigned int DF=0; DF<R.nDF(); DF++) {
		R.setInteractionDF(DF,1.0);
		mvec v = R.getReactionTo(&R);
		assert(v.size() == R.nDF()/R.nPhi);
		for(unsigned int i=0; i<v.size(); i++)
			the_ixn(i, DF/R.nPhi)[DF%R.nPhi] = v[i];
		pb.update(DF);
	}
	the_GF = makeIdentity(R.nDF()/R.nPhi, R.nPhi) - the_ixn;
}

void SymmetricSolver::calculateResult(ReactiveSet& R) {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	R.prepareIncident();
	R.finalState = R.incidentState;
	circulantMul(the_GF, R.finalState, R.nPhi);
	R.setFinalState(R.finalState);
	if(verbose) printf(" Done.\n");
}

void SymmetricSolver::selfInteract(ReactiveSet& R) {
	VarVec<double> v = R.incidentState;
	circulantMul(the_ixn, v, R.nPhi);
	R.setFinalState(v);
}

void SymmetricSolver::writeToFile(std::ostream& o) const {
	assert(false);
	//the_GF.writeToFile(o);
}

void SymmetricSolver::readFromFile(std::istream& s) const {
	assert(false);
	//the_GF.readFromFile(s);
}

void SymmetricSolver::cachedSolve(ReactiveSet& R, const std::string& fname) {
	
	if(!fname.size()) {
		solve(R);
		return;
	}
	
	std::ifstream ifs(fname.c_str(), std::ifstream::in | std::ifstream::binary);
	if(ifs.good()) {
		std::cout << "Loading previous SymmetricSolver solution from '" << fname << "'." << std::endl;
		readFromFile(ifs);
		ifs.close();
	} else {
		solve(R);
		std::ofstream ofs(fname.c_str(), std::ifstream::out | std::ifstream::binary);
		if(!ofs.good()) {
			std::cout << "Warning: SymmetricSolver unable to write to file '" << fname << "'." << std::endl;
			return;
		}
		writeToFile(ofs);
		ofs.close();
	}
}
