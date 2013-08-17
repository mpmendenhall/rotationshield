#include "GenericSolver.hh"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include <cassert>

GenericSolver::GenericSolver(): InteractionSolver(), the_GF(NULL) {}
GenericSolver::~GenericSolver() {
	if(the_GF) gsl_matrix_free(the_GF);
	the_GF = NULL;
}

void GenericSolver::solve(bool checkinversion) {
	
	buildInteractionMatrix();
	
	// invert matrix
	int sig;
	gsl_permutation* P = gsl_permutation_alloc(ndf);
	gsl_matrix* I = gsl_matrix_alloc(ndf,ndf);
	assert(!gsl_linalg_LU_decomp (the_GF, P, &sig));
	assert(!gsl_linalg_LU_invert (the_GF, P, I));
	gsl_permutation_free(P);
	gsl_matrix_free(the_GF);
	the_GF = I;

	//if(checkinversion)
	//	printf(" %g inversion error.",(double)BlockCMat<mdouble>::checkInversion(bc, the_GF));
}

void GenericSolver::buildInteractionMatrix() {
	the_GF = gsl_matrix_alloc(ndf,ndf);
	for(unsigned int i1=0; i1<N(); i1++) {
		for(unsigned int i2=0; i2<N(); i2++) {
			mmat mInt = getInteraction(i1,i2);
			for(unsigned int d1=0; d1<mInt.nRows(); d1++)
				for(unsigned int d2=0; d2<mInt.nCols(); d2++)
					gsl_matrix_set(the_GF,index(i1,d1),index(i2,d2), mInt(d1,d2));
		}
	}
}

void GenericSolver::calculateResult() {
	assert(the_GF);
	gsl_vector* inc = gsl_vector_alloc(ndf);
	gsl_vector* fin = gsl_vector_alloc(ndf);
	for(unsigned int i=0; i<ndf; i++)
		gsl_vector_set(inc,i,incidentState[i]);
	assert(!gsl_blas_dgemv(CblasNoTrans, 1., the_GF, inc, 0., fin));
	finalState = VarVec<mdouble>(ndf);
	for(unsigned int i=0; i<ndf; i++)
		finalState[i] = gsl_vector_get(fin,i);
	gsl_vector_free(inc);
	gsl_vector_free(fin);
	setSurfacelsOutputState();
}