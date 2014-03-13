#include "SymmetricSolver.hh"
#include "ProgressBar.hh"
#include <cassert>

SymmetricSolver::SymmetricSolver(): InteractionSolver() {
	verbose = true;
}

void SymmetricSolver::solve(ReactiveSet& R) {
	buildInteractionMatrix(R);
	printf("Solving for Green's Function... ");
	the_GF.invert();
	//if(checkinversion)
	//	printf(" %g inversion error.",(double)BlockCMat<mdouble>::checkInversion(bc, the_GF));
	printf(" Done.\n");
}

void SymmetricSolver::buildInteractionMatrix(ReactiveSet& R) {
	if(verbose) printf("Building interaction matrix for %i DF in groups of %i...\n", R.nDF(), R.nPhi);
	ProgressBar pb = ProgressBar(R.nDF(),R.nDF()/R.nPhi,verbose);
	the_GF = BlockCMat<mdouble>(R.nDF()/R.nPhi,R.nDF()/R.nPhi,R.nPhi);
	
	R.startInteractionScan();
	for(unsigned int DF=0; DF<R.nDF(); DF++) {
		R.setInteractionDF(DF);
		mvec v = R.getReactionTo(&R);
		assert(v.size() == R.nDF()/R.nPhi);
		for(unsigned int i=0; i<v.size(); i++)
			the_GF.getBlock(i, DF/R.nPhi)[DF%R.nPhi] = i*R.nPhi==DF ? 1-v[i] : -v[i];
		pb.update(DF);
	}
}

void SymmetricSolver::calculateResult(ReactiveSet& R) {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	R.setFinalState(the_GF * R.incidentState);
	if(verbose) printf(" Done.\n");
}
