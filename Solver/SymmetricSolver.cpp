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
	ProgressBar pb = ProgressBar(R.nDF()*R.nDF()/R.nPhi,R.nDF(),verbose);
	the_GF = BlockCMat<mdouble>(R.nDF()/R.nPhi,R.nDF()/R.nPhi,R.nPhi);
	
	unsigned int i;
	unsigned int j;
	R.startInteractionScan();
	for(unsigned int n=0; n<R.nDF()*R.nDF()/R.nPhi; n++) {
		mdouble v = R.nextInteractionTerm(i,j);
		the_GF.getBlock(i/R.nPhi, j/R.nPhi)[j%R.nPhi] = i==j ? 1-v : -v;
		pb.update(n);
	}
}

void SymmetricSolver::calculateResult(ReactiveSet& R) {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	R.setFinalState(the_GF * R.incidentState);
	if(verbose) printf(" Done.\n");
}
