#include "GenericSolver.hh"

void GenericSolver::solve(bool checkinversion) {
	
	buildInteractionMatrix();
	int errn;
	the_GF.invert(errn);
	//if(checkinversion)
	//	printf(" %g inversion error.",(double)BlockCMat<mdouble>::checkInversion(bc, the_GF));
}

void GenericSolver::buildInteractionMatrix() {
	the_GF = CLHEP::HepMatrix(ndf,ndf);
	for(unsigned int i1=0; i1<N(); i1++) {
		for(unsigned int i2=0; i2<N(); i2++) {
			mmat mInt = getInteraction(i1,i2);
			for(unsigned int d1=0; d1<mInt.nRows(); d1++)
				for(unsigned int d2=0; d2<mInt.nCols(); d2++)
					the_GF[index(i1,d1)][index(i2,d2)] = mInt(d1,d2);
		}
	}
}

void GenericSolver::calculateResult() {
	CLHEP::HepVector inc = CLHEP::HepVector(ndf);
	for(unsigned int i=0; i<ndf; i++)
		inc[i] = incidentState[i];
	CLHEP::HepVector fin = the_GF*inc;
	finalState = VarVec<mdouble>(ndf);
	for(unsigned int i=0; i<ndf; i++)
		finalState[i] = fin[i];
	setSurfacelsOutputState();
}