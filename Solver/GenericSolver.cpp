/* 
 * GenericSolver.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "GenericSolver.hh"
#include "ProgressBar.hh"

#include "gsl/gsl_vector.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include <cassert>

#ifdef WITH_LAPACKE
GenericSolver::GenericSolver(): InteractionSolver(true), my_SVD(NULL) {}

GenericSolver::~GenericSolver() {
	if(my_SVD) delete my_SVD;
}
#else
GenericSolver::GenericSolver(): InteractionSolver(true), the_GF(NULL) {}

GenericSolver::~GenericSolver() {
	if(the_GF) gsl_matrix_free(the_GF);
	the_GF = NULL;
}
#endif

void GenericSolver::solve(ReactiveSet& R) {
	
	buildInteractionMatrix(R);

#ifdef WITH_LAPACKE
	// calculate SVD of interaction matrix
	if(my_SVD) delete my_SVD;
	mmat GFI = mmat::identity(the_ixn.nRows()) - the_ixn;
	my_SVD = new LAPACKE_Matrix_SVD<double,double>(GFI);
	print_singular_values();
#else
	// invert matrix
	int sig;
	gsl_permutation* P = gsl_permutation_alloc(R.nDF());
	gsl_matrix* GFI = gsl_matrix_alloc(R.nDF(),R.nDF());
	assert(!gsl_linalg_LU_decomp (the_GF, P, &sig));
	assert(!gsl_linalg_LU_invert (the_GF, P, GFI));
	gsl_permutation_free(P);
	gsl_matrix_free(the_GF);
	the_GF = GFI;
#endif

}

void GenericSolver::buildInteractionMatrix(ReactiveSet& R) {
	if(verbose) printf("Building interaction matrix for %i DF...\n", R.nDF());
	ProgressBar pb = ProgressBar(R.nDF(),1,verbose);

#ifdef WITH_LAPACKE
	the_ixn = mmat(R.nDF(),R.nDF());
#else
	the_GF = gsl_matrix_alloc(R.nDF(),R.nDF());
#endif
	R.startInteractionScan();
	for(unsigned int DF=0; DF<R.nDF(); DF++) {
		R.setInteractionDF(DF,1.0);
		for(unsigned int phi=0; phi<R.nPhi; phi++) {
			mvec v = R.getReactionTo(&R,phi);
			assert(v.size() == R.nDF()/R.nPhi);
			for(unsigned int i=0; i<v.size(); i++)
#ifdef WITH_LAPACKE
				the_ixn(i*R.nPhi+phi, DF) = v[i];
#else
				gsl_matrix_set(the_GF, i*R.nPhi+phi, DF, (i*R.nPhi+phi==DF) ? 1-v[i] : -v[i]);
#endif
		}
		pb.update(DF);
	}
	R.setInteractionDF(R.nDF(),0);
}

void GenericSolver::calculateResult(ReactiveSet& R) {
#ifdef WITH_LAPACKE
	assert(my_SVD);
	mvec finalState = my_SVD->calc_pseudo_inverse(singular_epsilon) * R.incidentState;
#else
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
#endif
	
	R.setFinalState(finalState);
}
