#include "Integrator.hh"
#include <cassert>

bool integratingParams::verbose = false;

//
// 1D integrator
//

double generalIntegratingFunction(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end())
		return (it->second)[p->axis];
	mvec v = p->f(x,p->fparams);
	p->n_dim = v.size();
	p->m[x] = v;
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

double generalIntegratingFunctionNoCache(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	mvec v = p->f(x,p->fparams);
	p->n_dim = v.size();
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

void Integrator::setup_singularities(mdouble a, mdouble b) {
	_singularities.clear();
	sng_0 = sng_1 = singularities.size();
	if(!singularities.size()) return;
	
	//printf("Checking %i singularities in %g--%g...\n",singularities.size(),a,b);
	std::vector<double>::iterator sa = std::lower_bound(singularities.begin(),singularities.end(),a);
	std::vector<double>::iterator sb = std::upper_bound(singularities.begin(),singularities.end(),b);
	sng_0 = sa-singularities.begin();
	sng_1 = sb-singularities.begin();
	//printf("Range = %i--%i\n",sng_0,sng_1);
	
	_singularities.push_back(a);
	if(sa != singularities.end() && *sa != a)
		_singularities.push_back(a);
	_singularities.insert(_singularities.end(),sa,sb);
	if(_singularities.back() != b)
		_singularities.push_back(b);
	//for(unsigned int i=0; i<_singularities.size(); i++)
	//	printf("\t%g\n",_singularities[i]);
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
	
	//return _integrate(p, &generalIntegratingFunction, a, b) * (isneg?-1:1);
	return _integrate(p, &generalIntegratingFunctionNoCache, a, b) * (isneg?-1:1);
}

mvec Integrator::_integrate(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b) {
	p.m.clear();

	gsl_function F;
	F.function = integf;
	F.params = &p;
	double r,e;
	int er;
	
	mvec v;
	p.axis = 0;
	do {
		if(adaptive) {
			if(_singularities.size())
				er = gsl_integration_qagp(&F, &_singularities[0], _singularities.size(), abs_err, rel_err, 1024, gslIntegrationWS, &r, &e);
			else
				er = gsl_integration_qags(&F, a, b, abs_err, rel_err, 1024, gslIntegrationWS, &r, &e);
			//er = gsl_integration_qag(&F, a, b, abs_err, rel_err, 1024, GSL_INTEG_GAUSS15, gslIntegrationWS, &r, &e);
			if(er) printf("(*Integration Warning %s %i\t(axis %i = %g:\t%g err)*)\n",p.nm.c_str(),er,p.axis,r,e);
		} else {
			size_t neval;
			er = gsl_integration_qng(&F, a, b, abs_err, rel_err, &r, &e, &neval);
			if(er) printf("(*Integration Warning %s %i\t(axis %i = %g:\t%i evals, %g err)*)\n",p.nm.c_str(),er,p.axis,r,(int)neval,e);
		}
		v.push_back(r);
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


/// Contains arguments for general 2D integrating function
class integratingParams_2D: public integratingParams {
public:
	mvec (*f2)(mdouble,mdouble,void *);	//< pointer to the 2D function being integrated
	Integrator yIntegrator;				//< Integrator for y-direction integrals
	mdouble x;							//< x value being integrated
	mdouble y0;							//< y lower bound of integration
	mdouble y1;							//< y upper bound of integration
};

// slice of a 2D function at constant x
mvec xslice(mdouble y, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	return p->f2(p->x,y,p->fparams);
}

double generalIntegratingFunction2D(double x, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) {
		//if(p->verbose) std::cout << x << " " << p->axis << " " << it->second << std::endl;
		return (it->second)[p->axis];
	}
	//std::cout << "x = " << x << std::endl;
	p->x = x;
	mvec v = p->yIntegrator.integrate(&xslice, p->y0, p->y1, params);
	p->n_dim = v.size();
	p->m[x] = v;
	//if(p->verbose) std::cout << x << v << std::endl;
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

void Integrator2D::setup_singularities(mdouble a, mdouble b) {
	singularities.clear();
	for(std::vector<vec2>::iterator it = xysingularities.begin(); it != xysingularities.end(); it++)
		singularities.push_back((*it)[0]);
	Integrator::setup_singularities(a,b);
}

mvec Integrator2D::integrate(mvec (*f)(mdouble,mdouble,void*), mdouble x0, mdouble x1, mdouble y0, mdouble y1, void* params) {
	
	integratingParams_2D p;
	p.f2 = f;
	p.fparams = params;
	p.nm = "2D Integral";

	setup_singularities(x0,x1);
	for(unsigned int i=sng_0; i<sng_1; i++)
		p.yIntegrator.singularities.push_back(xysingularities[i][1]);
		
	p.y0 = y0;
	p.y1 = y1;
	return _integrate(p, &generalIntegratingFunction2D, x0, x1);
}

