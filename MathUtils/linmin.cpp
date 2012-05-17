#include "linmin.hh"

CLHEP::HepVector lsmin(CLHEP::HepMatrix coeffs, CLHEP::HepVector rslt) {
	CLHEP::HepSymMatrix M = CLHEP::HepDiagMatrix(coeffs.num_row(),1).similarityT(coeffs);
	CLHEP::HepVector v = coeffs.T() * rslt;
	return solve(M,v);
}

double ssresids(CLHEP::HepMatrix coeffs, CLHEP::HepVector x, CLHEP::HepVector y) {
	return (coeffs * x - y).normsq();
}

CLHEP::HepVector polyQuadraticFit(CLHEP::HepMatrix coords, CLHEP::HepVector values, double& resids) {
	int nvars = coords.num_col();
	int nparams = 1 + nvars + (nvars+1)*nvars/2; //# of parameters
	CLHEP::HepMatrix coeffs(coords.num_row(),nparams);
	// build coefficients matrix
	for(int i=0; i<coords.num_row(); i++)
	{
		coeffs[i][0] = 1.0;
		for(int j=0; j<nvars; j++)
			coeffs[i][j+1] = coords[i][j];
		int n=0;
		for(int j=0; j<nvars; j++)
			for(int k=0; k<=j; k++)
				coeffs[i][1+nvars+(n++)] = coords[i][j]*coords[i][k];
	}
	CLHEP::HepVector fitv = lsmin(coeffs,values);
	resids = sqrt(ssresids(coeffs,fitv,values)/mdouble(coords.num_row()));
	return fitv;
}


double polynomialFit(CLHEP::HepMatrix coords, CLHEP::HepVector values, Polynomial<3,mdouble>& p) {
	int nparams = p.terms.size();
	
	// build coefficients matrix
	CLHEP::HepMatrix coeffs(coords.num_row(),nparams);
	Vec<3,mdouble> coord;
	p.it = p.terms.begin();
	for(int j=0; j<nparams; j++) {
		Monomial<3,mdouble,unsigned int> m = Monomial<3,mdouble,unsigned int>(1.0,p.it->first);
		for(int i=0; i<coords.num_row(); i++)
		{
			for(int c=0; c<3; c++) coord[c] = coords[i][c];
			coeffs[i][j] = m(coord);
		}
		p.it++;
	}
	
	// solve, store coeffs to polynomial, and return
	CLHEP::HepVector fitv = lsmin(coeffs,values);
	p.it = p.terms.begin();
	for(int j=0; j<nparams; j++) {
		p.terms[p.it->first] = fitv[j];
		p.it++;
	}
	return sqrt(ssresids(coeffs,fitv,values)/mdouble(coords.num_row()));
}


void test_fitting() {
	CLHEP::HepMatrix coords(9*9*9,3);
	CLHEP::HepVector v(9*9*9);
	double x,y,z;
	double r;
	for(int i=0; i<9*9*9; i++)
	{
		x = i/81; y = (i%81)/9; z=i%9;
		coords[i][0] = x;
		coords[i][1] = y;
		coords[i][2] = z;
		v[i] = 3.14 + 1*x + 2*y + 3*z +4*x*x + 5*y*x + 6*y*y + 7*z*x + 8*z*y + 9*z*z;
	}
	std::cout << coords << v << polyQuadraticFit(coords,v,r);
}

