#include "linmin.hh"
#include <cassert>
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"


gsl_vector* lsmin(gsl_matrix* coeffs, const gsl_vector* rslt, gsl_vector* resid) {
	
	assert(coeffs && rslt && resid);
	assert(coeffs->size1 >= coeffs->size2);
	assert(resid->size == coeffs->size1);
	assert(rslt->size == coeffs->size1);
	
	gsl_vector* tau = gsl_vector_calloc(coeffs->size2);
	assert(!gsl_linalg_QR_decomp(coeffs,tau));
	gsl_vector* x = gsl_vector_calloc(coeffs->size2);
	assert(!gsl_linalg_QR_lssolve(coeffs, tau, rslt, x, resid));
	
	gsl_vector_free(tau);
	gsl_matrix_free(coeffs);
	return x;
}

double polynomialFit(const gsl_matrix* coords, const gsl_vector* values, Polynomial<3,double>& p) {
	int nparams = p.terms.size();
	assert(coords && values);
	assert((unsigned int)nparams <= values->size);
	assert(coords->size1 == values->size);
	assert(coords->size2 == 3);
	
	// build coefficients matrix
	gsl_matrix* coeffs = gsl_matrix_alloc(coords->size1,nparams);
	Vec<3,double> coord;
	p.it = p.terms.begin();
	for(int j=0; j<nparams; j++) {
		Monomial<3,double,unsigned int> m = Monomial<3,double,unsigned int>(1.0,p.it->first);
		for(unsigned int i=0; i<values->size; i++) {
			for(int c=0; c<3; c++) coord[c] = gsl_matrix_get(coords,i,c);
			gsl_matrix_set(coeffs,i,j,m(coord));
		}
		p.it++;
	}
	
	// fit, cleanup, return
	gsl_vector* resid = gsl_vector_calloc(values->size);
	gsl_vector* fitv = lsmin(coeffs,values,resid);
	p.it = p.terms.begin();
	for(int j=0; j<nparams; j++) {
		p.terms[p.it->first] = gsl_vector_get(fitv,j);
		p.it++;
	}
	double rsresid =  gsl_blas_dnrm2(resid);
	gsl_vector_free(fitv);
	gsl_vector_free(resid);
	return rsresid/sqrt(values->size);
}
