#include "SymmetricSolver.hh"
#include "ProgressBar.hh"
#include "PathUtils.hh"
#include <cassert>
#include <fstream>

void SymmetricSolver::circulantMul(const BlockCMat& M, mvec& v, unsigned int nPhi) {
	assert(!(v.size()%nPhi));
	assert(M.nCols()*nPhi == v.size());
	
	// stuff vector into vector-of-vectors for circulant blocks
	VarVec<mvec> vv;
	for(unsigned int i=0; i<v.size()/nPhi; i++)
		vv.push_back(mvec(&v[i*nPhi], &v[i*nPhi]+nPhi));
	vv = M.lMultiply<mvec,mvec>(vv);
	
	// pull data back out
	v.getData().resize(M.nRows()*nPhi);
	for(unsigned int i=0; i<M.nRows(); i++)
		for(unsigned int j=0; j<nPhi; j++)
			v[i*nPhi+j] = vv[i][j];
}

double SymmetricSolver::checkInversion(const BlockCMat& M, const BlockCMat& MI, unsigned int nPhi) {
	return (M * MI - makeBlockCMatIdentity(M.nRows(),nPhi)).getData().max_norm_L2();
}

void SymmetricSolver::solve(ReactiveSet& R) {
	buildInteractionMatrix(R);
	printf("Solving for Green's Function... ");
	if(the_GF) delete(the_GF);
	BlockCMat ImR = makeBlockCMatIdentity(R.nDF()/R.nPhi, R.nPhi) - the_ixn;
	the_GF = new BlockCMat_SVD(ImR);
	std::cout << "Inversion error " << checkInversion(ImR, the_GF->calc_pseudo_inverse(singular_epsilon), R.nPhi) << "\n";
	print_singular_values();
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
}

void SymmetricSolver::calculateResult(ReactiveSet& R) {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	R.prepareIncident();
	R.finalState = R.incidentState;
	circulantMul(the_GF->calc_pseudo_inverse(singular_epsilon), R.finalState, R.nPhi);
	R.setFinalState(R.finalState);
	if(verbose) printf(" Done.\n");
}

void SymmetricSolver::selfInteract(ReactiveSet& R) {
	R.finalState = R.incidentState;
	circulantMul(the_ixn, R.finalState, R.nPhi);
	R.setFinalState(R.finalState);
}

void SymmetricSolver::print_singular_values() const {
	if(!the_GF) return;
	std::cout << "Most Singular values:\n";
	const std::vector<double>& S = the_GF->singular_values();
	for(unsigned int i=0; i<S.size(); i++) {
		std::cout << S[i]/S.back() << "\t";
		if(S[i] > 0.01*S.back()) break;
	}
	std::cout << "\n\n";
}

VarVec<double> SymmetricSolver::get_singular_vector(unsigned int i) const {
	if(the_GF) return the_GF->getRightSVec(i);
	return VarVec<double>();
}

double SymmetricSolver::get_singular_value(unsigned int i) const {
	if(!the_GF) return 0;
	const std::vector<double>& S = the_GF->singular_values();
	assert(i<S.size());
	return S[i]/S.back();
}

void SymmetricSolver::writeToFile(std::ostream& o) const {
	writeString("(SymmetricSolver)",o);
	the_ixn.writeToFile(o);
	o.write((char*)&singular_epsilon,		sizeof(singular_epsilon));
	o.write((char*)&the_GF,					sizeof(the_GF));
	if(the_GF) the_GF->writeToFile(o);
	writeString("(/SymmetricSolver)",o);
}

SymmetricSolver* SymmetricSolver::readFromFile(std::istream& s) {
	SymmetricSolver* foo = new SymmetricSolver();
	checkString("(SymmetricSolver)",s);
	foo->the_ixn = BlockCMat::readFromFile(s);
	s.read((char*)&foo->singular_epsilon,	sizeof(foo->singular_epsilon));
	s.read((char*)&foo->the_GF,				sizeof(foo->the_GF));
	if(foo->the_GF) foo->the_GF = BlockCMat_SVD::readFromFile(s);
	checkString("(/SymmetricSolver)",s);
	return foo;
}

SymmetricSolver* SymmetricSolver::cachedSolve(ReactiveSet& R, const std::string& fname) {
	
	SymmetricSolver* foo = NULL;
	if(!fname.size()) {
		foo = new SymmetricSolver();
		foo->solve(R);
		return foo;
	}
	
	std::string ext_fname = fname + "_" + std::to_string(R.nDF()/R.nPhi) + "_" + std::to_string(R.nPhi) + ".slvdat";
	makePath(ext_fname, true);
	
	std::ifstream ifs(ext_fname.c_str(), std::ifstream::in | std::ifstream::binary);
	if(ifs.good()) {
		std::cout << "Loading previous SymmetricSolver solution from '" << ext_fname << "'." << std::endl;
		foo = SymmetricSolver::readFromFile(ifs);
		ifs.close();
	} else {
		std::cout << "No previous solution found at '" << ext_fname << "'; Solving." << std::endl;
		foo = new SymmetricSolver();
		foo->solve(R);
		std::cout << "Saving SymmetricSolver solution to '" << ext_fname << "'." << std::endl;
		std::ofstream ofs(ext_fname.c_str(), std::ofstream::out | std::ofstream::binary);
		if(!ofs.good()) {
			std::cout << "Warning: SymmetricSolver unable to write to file '" << ext_fname << "'." << std::endl;
			return foo;
		}
		foo->writeToFile(ofs);
		ofs.close();
	}
	
	return foo;
}
