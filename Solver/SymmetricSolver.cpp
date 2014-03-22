#include "SymmetricSolver.hh"
#include "ProgressBar.hh"
#include <cassert>
#include <fstream>

void SymmetricSolver::solve(ReactiveSet& R) {
	buildInteractionMatrix(R);
	printf("Solving for Green's Function... ");
	the_GF.invert();
	//if(checkinversion)
	//	printf(" %g inversion error.",(double)BlockCMat<mdouble>::checkInversion(bc, the_GF));
	printf(" Done.\n");
}

void SymmetricSolver::buildInteractionMatrix(ReactiveSet& R) {
	if(verbose) printf("Building interaction matrix for %i = %i x %i DF...\n", R.nDF(), R.nPhi, R.nDF()/R.nPhi);
	ProgressBar pb = ProgressBar(R.nDF(), R.nPhi, verbose);
	the_ixn = BlockCMat<mdouble>(R.nDF()/R.nPhi,R.nDF()/R.nPhi,R.nPhi);
	
	R.startInteractionScan();
	for(unsigned int DF=0; DF<R.nDF(); DF++) {
		R.setInteractionDF(DF,1.0);
		mvec v = R.getReactionTo(&R);
		assert(v.size() == R.nDF()/R.nPhi);
		for(unsigned int i=0; i<v.size(); i++)
			the_ixn.getBlock(i, DF/R.nPhi)[DF%R.nPhi] = v[i];
		pb.update(DF);
	}
	the_GF = BlockCMat<mdouble>::identity(R.nDF()/R.nPhi,R.nPhi) - the_ixn;
}

void SymmetricSolver::calculateResult(ReactiveSet& R) {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	R.prepareIncident();
	R.setFinalState(the_GF * R.incidentState);
	if(verbose) printf(" Done.\n");
}

void SymmetricSolver::selfInteract(ReactiveSet& R) {
	R.setFinalState(the_ixn * R.incidentState);
}

void SymmetricSolver::writeToFile(std::ostream& o) const {
	the_GF.writeToFile(o);
}

void SymmetricSolver::readFromFile(std::istream& s) const {
	the_GF.readFromFile(s);
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
