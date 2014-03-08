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

mvec Integrator::integrate(mvec (*f)(mdouble,void*), mdouble a, mdouble b, void* params) {
	integratingParams p;
	p.fparams = params;
	p.f = f;
	return _integrate(p, &generalIntegratingFunction, a, b);
}

mvec Integrator::_integrate(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b) {
	p.m.clear();

	gsl_function F;
	F.function = integf;
	F.params = &p;
	double r,e;
	size_t neval;
	
	mvec v;
	p.axis = 0;
	do {
		int er = gsl_integration_qng(&F, a, b, 1e-7, 1e-7, &r, &e, &neval);
		if(er)
			printf("(*INTEGRATION WARNING*)\n");
		v.push_back(r);
		p.axis++;
	} while (p.axis < p.n_dim);
	
	return v;
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
	p->x = x;
	mvec v = p->yIntegrator.integrate(&xslice, p->y0, p->y1, params);
	p->n_dim = v.size();
	p->m[x] = v;
	//if(p->verbose) std::cout << x << v << std::endl;
	assert(p->axis < p->n_dim);
	return v[p->axis];
}

mvec Integrator2D::integrate(mvec (*f)(mdouble,mdouble,void*), mdouble x0, mdouble x1, mdouble y0, mdouble y1, void* params) {
	integratingParams_2D p;
	p.f2 = f;
	p.fparams = params;
	p.y0 = y0;
	p.y1 = y1;
	
	return _integrate(p, &generalIntegratingFunction2D, x0, x1);
}

