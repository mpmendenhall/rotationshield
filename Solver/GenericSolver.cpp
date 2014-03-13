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
	the_GF = gsl_matrix_alloc(R.nDF(),R.nDF());
	R.startInteractionScan();
	for(unsigned int DF=0; DF<R.nDF(); DF++) {
		R.setInteractionDF(DF,1.0);
		for(unsigned int phi=0; phi<R.nPhi; phi++) {
			mvec v = R.getReactionTo(&R,phi);
			assert(v.size() == R.nDF()/R.nPhi);
			for(unsigned int i=0; i<v.size(); i++)
				gsl_matrix_set(the_GF, i*R.nPhi+phi, DF, (i*R.nPhi+phi==DF) ? 1-v[i] : -v[i]);
		}
	}
}

void GenericSolver::calculateResult(ReactiveSet& R) {
	R.prepareIncident();
	
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