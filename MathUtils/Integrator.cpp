#include "Integrator.hh"
#include <cassert>

bool integratingParams::verbose = false;

//
// 1D integrator
//

void Integrator::setup_singularities(mdouble a, mdouble b) {
	_singularities.clear();
	if(!singularities.size()) return;
	
	std::set<double>::iterator sa = std::lower_bound(singularities.begin(),singularities.end(),a);
	std::set<double>::iterator sb = std::upper_bound(singularities.begin(),singularities.end(),b);
	
	_singularities.push_back(a);
	if(sa != singularities.end() && *sa != a)
		_singularities.push_back(a);
	_singularities.insert(_singularities.end(),sa,sb);
	if(_singularities.back() != b)
		_singularities.push_back(b);
}

double Integrator::integrate(double (*f)(double, void*), double a, double b, void* params) {
	
	bool isneg = b<a;
	if(isneg) std::swap(a,b);
	setup_singularities(a,b);
	
	gsl_function F;
	F.function = f;
	F.params = params;
	
	return _integrate(&F,a,b)*(isneg?-1:1);
}

double Integrator::_integrate(gsl_function* F, double a, double b) {

	double r=0;
	double e;
	int er = 0;
	size_t neval;
		
	if(myMethod == INTEG_GSL_QNG) {
		er = gsl_integration_qng(F, a, b, abs_err, rel_err, &r, &e, &neval);
		if(er) printf("(*Integration Warning %i\t(%g:\t%i evals, %g err)*)\n",er,r,(int)neval,e);
	} else {
		if(myMethod == INTEG_GSL_QAG)
			er = gsl_integration_qag(F, a, b, abs_err, rel_err, 256, GSL_INTEG_GAUSS15, gslIntegrationWS, &r, &e);
		else if (myMethod == INTEG_GSL_CQUAD)
			er = gsl_integration_cquad(F, a, b, abs_err, rel_err, gsl_cqd_ws, &r, &e, &neval);
		else if(myMethod == INTEG_GSL_QAGS)
			er = gsl_integration_qags(F, a, b, abs_err, rel_err, 256, gslIntegrationWS, &r, &e);
		else if(myMethod == INTEG_GSL_QAGP) {
			if(_singularities.size() < 2) {
				_singularities.clear();
				_singularities.push_back(a);
				_singularities.push_back(b);
			}
			er = gsl_integration_qagp(F, &_singularities[0], _singularities.size(), abs_err, rel_err, 256, gslIntegrationWS, &r, &e);
		} else { assert(false); }
		if(er) printf("(*Integration Warning %i\t(%g:\t%g err)*)\n",er,r,e);
	}
	return r;
}

//
// 1D multivariate, cacheing point results on other axes
//

double generalIntegratingFunction(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) return (it->second)[p->axis];
	
	mvec v = p->f(x,p->fparams);
	p->n_dim = v.size();
	p->m[x] = v;
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

double generalIntegratingFunctionNoCache(double x, void* params) {
	// non-cacheing version always re-evaluates function components
	integratingParams* p = (integratingParams*)params;
	mvec v = p->f(x,p->fparams);
	p->n_dim = v.size();
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

mvec Integrator::integrate(mvec (*f)(mdouble,void*), mdouble a, mdouble b, void* params) {

	integratingParams p;
	p.fparams = params;
	p.f = f;
	p.nm = "1D Integral";
	
	bool isneg = b<a;
	if(isneg) std::swap(a,b);
	
	// set up singularities list
	setup_singularities(a,b);
	
	return _integrate_v(p, &generalIntegratingFunction, a, b) * (isneg?-1:1);
	//return _integrate_v(p, &generalIntegratingFunctionNoCache, a, b) * (isneg?-1:1);
}

mvec Integrator::_integrate_v(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b) {
	p.m.clear();
	
	gsl_function F;
	F.function = integf;
	F.params = &p;
		
	mvec v;
	p.axis = 0;
	do {
		v.push_back(_integrate(&F,a,b));
		p.axis++;
	} while (p.axis < p.n_dim);
	
	return v;
}

void Integrator::printErrorCodes() {
	printf("GSL Integrator error codes:\n");
	printf("%i: the maximum number of subdivisions was exceeded.\n",GSL_EMAXITER);
    printf("%i: cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.\n",GSL_EROUND);
    printf("%i: a non-integrable singularity or other bad integrand behavior was found in the integration interval.\n",GSL_ESING);
	printf("%i: the integral is divergent, or too slowly convergent to be integrated numerically.\n",GSL_EDIVERGE);
}





//
// 2D integrator
//


void Integrator2D::setup_singularities(mdouble a, mdouble b) {
	singularities.clear();
	for(std::vector<vec2>::iterator it = xysingularities.begin(); it != xysingularities.end(); it++)
		singularities.insert((*it)[0]);
	Integrator::setup_singularities(a,b);
}

/// Contains arguments for general 2D integrating function
class integratingParams_2D {
public:
	mdouble (*f2)(vec2,void*);	//< pointer to the 2D function being integrated
	void* fparams;				//< parameters for integrating function
	Integrator* yIntegrator;	//< Integrator for y-direction integrals
	mdouble x;					//< x value being integrated
	mdouble y0;					//< y lower bound of integration
	mdouble y1;					//< y upper bound of integration
};

// slice of a 2D function at constant x
double xslice(double y, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	return p->f2(vec2(p->x,y),p->fparams);
}
// integral across slice of a 2D function
double xslice_integral(double x, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	p->x = x;
	return p->yIntegrator->integrate(&xslice, p->y0, p->y1, p);
}

mdouble Integrator2D::integrate2D(mdouble (*f)(vec2,void*), vec2 ll, vec2 ur, void* params) {
		
	integratingParams_2D p;
	p.yIntegrator = &yIntegrator;
	p.f2 = f;
	p.fparams = params;
	p.y0 = ll[1];
	p.y1 = ur[1];
	
	setup_singularities(ll[0],ur[0]); // TODO better selection of y singularities
	yIntegrator.singularities.clear();
	for(std::vector<vec2>::iterator it = xysingularities.begin(); it != xysingularities.end(); it++)
		yIntegrator.singularities.insert((*it)[1]);
	
	return Integrator::integrate(&xslice_integral, ll[0], ur[0], &p);
}

// multi-variate

/// Contains arguments for general 2D integrating function
class integratingParams_V2D: public integratingParams {
public:
	mvec (*f2)(vec2,void *);	//< pointer to the 2D function being integrated
	Integrator* yIntegrator;	//< Integrator for y-direction integrals
	mdouble x;					//< x value being integrated
	mdouble y0;					//< y lower bound of integration
	mdouble y1;					//< y upper bound of integration
};

// slice of a 2D function at constant x
mvec xslice_v(mdouble y, void* params) {
	integratingParams_V2D* p = (integratingParams_V2D*)params;
	return p->f2(vec2(p->x,y),p->fparams);
}


double generalIntegratingFunction2D(double x, void* params) {

	integratingParams_V2D* p = (integratingParams_V2D*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) { return (it->second)[p->axis]; }

	integratingParams_V2D p2;
	p2.f2 = p->f2;
	p2.fparams = p->fparams;
	p2.x = x;
	p2.yIntegrator = p->yIntegrator;
	p2.y0 = p->y0;
	p2.y1 = p->y1;
	mvec v = p->yIntegrator->integrate(&xslice_v, p->y0, p->y1, &p2);
	
	p->n_dim = v.size();
	p->m[x] = v;
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

mvec Integrator2D::integrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* params) {
	
	integratingParams_V2D p;
	p.yIntegrator = &yIntegrator;
	p.f2 = f;
	p.fparams = params;
	p.nm = "2D Integral";
	
	setup_singularities(ll[0],ur[0]); // TODO better selection of y singularities
	yIntegrator.singularities.clear();
	for(std::vector<vec2>::iterator it = xysingularities.begin(); it != xysingularities.end(); it++)
		yIntegrator.singularities.insert((*it)[1]);
	
	p.y0 = ll[1];
	p.y1 = ur[1];
	return _integrate_v(p, &generalIntegratingFunction2D, ll[0], ur[0]);
}

