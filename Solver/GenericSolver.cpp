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

void GenericSolver::solve(ReactiveSet& R) {
	
	buildInteractionMatrix(R);
	
	// invert matrix
	int sig;
	gsl_permutation* P = gsl_permutation_alloc(R.nDF());
	gsl_matrix* I = gsl_matrix_alloc(R.nDF(),R.nDF());
	assert(!gsl_linalg_LU_decomp (the_GF, P, &sig));
	assert(!gsl_linalg_LU_invert (the_GF, P, I));
	gsl_permutation_free(P);
	gsl_matrix_free(the_GF);
	the_GF = I;
}

void GenericSolver::buildInteractionMatrix(ReactiveSet& R) {
	/*
	the_GF = gsl_matrix_alloc(R.nDF(),R.nDF());
	unsigned int i;
	unsigned int j;
	R.startInteractionScan();
	for(unsigned int n=0; n<R.nDF()*R.nDF(); n++) {
		mdouble v = R.nextInteractionTerm(i,j);
		gsl_matrix_set(the_GF, i, j, i==j ? 1-v : -v);
	}
	*/
}

void GenericSolver::calculateResult(ReactiveSet& R) {
	assert(the_GF);
	gsl_vector* inc = gsl_vector_alloc(R.nDF());
	gsl_vector* fin = gsl_vector_alloc(R.nDF());
	for(unsigned int i=0; i<R.nDF(); i++)
		gsl_vector_set(inc,i,R.incidentState[i]);
	assert(!gsl_blas_dgemv(CblasNoTrans, 1., the_GF, inc, 0., fin));
	mvec finalState = mvec(R.nDF());
	for(unsigned int i=0; i<R.nDF(); i++)
		finalState[i] = gsl_vector_get(fin,i);
	gsl_vector_free(inc);
	gsl_vector_free(fin);
	
	R.setFinalState(finalState);
}