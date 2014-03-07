#include "Integrator.hh"


bool integratingParams::verbose = false;

//
// 1D integrator
//

double generalIntegratingFunction(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	std::map<double,vec3>::iterator it = p->m.find(x);
	if(it != p->m.end()) { if(p->verbose) std::cout << x << " " << p->axis << " " << it->second << std::endl; return (it->second)[p->axis]; }
	vec3 v = p->f(x,p->fparams);
	p->m[x] = v;
	if(p->verbose) std::cout << x << v << std::endl;
	return v[p->axis];
}

vec3 Integrator::integrate(vec3 (*f)(mdouble,void*),mdouble a, mdouble b, void* params) {
	integratingParams p;
	p.fparams = params;
	p.f = f;
	return _integrate(p, &generalIntegratingFunction, a, b);
}

vec3 Integrator::_integrate(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b) {
	p.m.clear();

	gsl_function F;
	F.function = integf;
	F.params = &p;
	double r,e;
	size_t neval;
	
	vec3 v;
	for(int i=0; i<3; i++) {
		p.axis = i;
		int er = gsl_integration_qng(&F, a, b, 1e-7, 1e-7, &r, &e, &neval);
		if(er)
			printf("(*INTEGRATION WARNING*)\n");
		v[i] = r;
	}
	return v;
}


//
// 2D integrator
//


/// Contains arguments for general 2D integrating function
class integratingParams_2D: public integratingParams {
public:
	vec3 (*f2)(mdouble,mdouble,void *);	//< pointer to the 2D function being integrated
	Integrator yIntegrator;				//< Integrator for y-direction integrals
	mdouble x;							//< x value being integrated
	mdouble y0;							//< y lower bound of integration
	mdouble y1;							//< y upper bound of integration
};

// slice of a 2D function at constant x
vec3 xslice(mdouble y, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	return p->f2(p->x,y,p->fparams);
}

double generalIntegratingFunction2D(double x, void* params) {
	integratingParams_2D* p = (integratingParams_2D*)params;
	std::map<double,vec3>::iterator it = p->m.find(x);
	if(it != p->m.end()) { if(p->verbose) std::cout << x << " " << p->axis << " " << it->second << std::endl; return (it->second)[p->axis]; }
	p->x = x;
	vec3 v = p->yIntegrator.integrate(&xslice, p->y0, p->y1, params);
	p->m[x] = v;
	if(p->verbose) std::cout << x << v << std::endl;
	return v[p->axis];
}

vec3 Integrator2D::integrate(vec3 (*f)(mdouble,mdouble,void*), mdouble x0, mdouble x1, mdouble y0, mdouble y1, void* params) {
	integratingParams_2D p;
	p.f2 = f;
	p.fparams = params;
	p.y0 = y0;
	p.y1 = y1;
	
	return _integrate(p, &generalIntegratingFunction2D, x0, x1);
}

