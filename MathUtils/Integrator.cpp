#include "Integrator.hh"

double generalIntegratingFunction(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	std::map<double,vec3>::iterator it = p->m.find(x);
	if(it != p->m.end()) { if(p->verbose) std::cout << x << " " << p->axis << " " << it->second << std::endl; return (it->second)[p->axis]; }
	vec3 v = p->f(x,p->fparams);
	p->m[x] = v;
	if(p->verbose) std::cout << x << v << std::endl;
	return v[p->axis];
}