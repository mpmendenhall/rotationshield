#include "SymmetricSolver.hh"
#include "ProgressBar.hh"

SymmetricSolver::SymmetricSolver(ShieldBuilder* g): InteractionSolver(g->nTheta), nZ(g->nSegments()), nTheta(g->nTheta) {
	verbose = true;
	for(unsigned int z=0; z<nZ; z++)
		for(unsigned int th=0; th<nTheta; th++)
			addSurfacel(g->genElement(z, th));
}

void SymmetricSolver::solve(bool checkinversion) {
	
	buildInteractionMatrix();
	printf("Solving for Green's Function... ");
	the_GF.invert();
	printf(" Done.\n");
}

void SymmetricSolver::buildInteractionMatrix() {
	if(verbose) printf("Building interaction matrix for %i DF in groups of %i...\n",ndf,groupSize);
	ProgressBar pb = ProgressBar(nZ,1,verbose);
	the_GF = BlockCMat<mdouble>(ndf/nTheta,ndf/nTheta,nTheta);
	for(unsigned int z0=0; z0<nZ; z0++) {
		pb.update(z0);
		for(unsigned int z1 = 0; z1 < nZ; z1++) {
			for(unsigned int th = 0; th<nTheta; th++) {
				mmat mInt = getInteraction(z0*nTheta,z1*nTheta+th);
				for(unsigned int d1=0; d1<mInt.nRows(); d1++)
					for(unsigned int d2=0; d2<mInt.nCols(); d2++)
						the_GF.getBlock(index(z0*nTheta,d1)/nTheta,index(z1*nTheta,d2)/nTheta)[th] = mInt(d1,d2);
			}
		}
	}
}

void SymmetricSolver::calculateResult() {
	if(verbose) printf("Calculating resulting surface current..."); fflush(stdout);
	finalState = the_GF * incidentState;
	setSurfacelsOutputState();
	if(verbose) printf(" Done.\n");
}
