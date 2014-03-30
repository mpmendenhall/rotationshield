/* 
 * Integrator.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * this code includes wrappers for integrators from the Gnu Scientific Library,
 * https://www.gnu.org/software/gsl/
 * and wrappers for multi-dimensional integration routines in
 * "cubature," (c) 2005-2013 Steven G. Johnson,
 * http://ab-initio.mit.edu/cubature/
 * re-distributed with in this project under GPL v2 or later;
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

#include "Integrator.hh"
#include <cassert>
#include "Angles.hh"

bool integratingParams::verbose = false;

//
// 1D integrator
//

Integrator::Integrator():
rel_err(1e-4), abs_err(1e-5), err_count(0), myMethod(INTEG_GSL_QNG),
gslIntegrationWS(gsl_integration_workspace_alloc(INTEG_WS_SIZE)),
gsl_cqd_ws(gsl_integration_cquad_workspace_alloc(INTEG_WS_SIZE/4)) {
	gsl_set_error_handler_off();
}

Integrator::~Integrator() {
	gsl_integration_workspace_free(gslIntegrationWS);
	gsl_integration_cquad_workspace_free(gsl_cqd_ws);
}


void Integrator::setup_singularities(double a, double b) {
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
		if(er) {
			printf("(*Integration Error %i '%s'\t(%g:\t%i evals, %g err)*)\n",er,gsl_strerror(er),r,(int)neval,e);
			err_count++;
		}
	} else {
		if(myMethod == INTEG_GSL_QAG)
			er = gsl_integration_qag(F, a, b, abs_err, rel_err, INTEG_WS_SIZE, GSL_INTEG_GAUSS15, gslIntegrationWS, &r, &e);
		else if (myMethod == INTEG_GSL_CQUAD)
			er = gsl_integration_cquad(F, a, b, abs_err, rel_err, gsl_cqd_ws, &r, &e, &neval);
		else if(myMethod == INTEG_GSL_QAGS)
			er = gsl_integration_qags(F, a, b, abs_err, rel_err, INTEG_WS_SIZE, gslIntegrationWS, &r, &e);
		else if(myMethod == INTEG_GSL_QAGP) {
			if(_singularities.size() < 2) {
				_singularities.clear();
				_singularities.push_back(a);
				_singularities.push_back(b);
			}
			er = gsl_integration_qagp(F, &_singularities[0], _singularities.size(), abs_err, rel_err, INTEG_WS_SIZE, gslIntegrationWS, &r, &e);
		} else { assert(false); }
		if(er) {
			printf("(*Integration Warning %i '%s'\t(%g:\t%g err)*)\n",er,gsl_strerror(er),r,e);
			err_count++;
		}
	}
	return r;
}

//
// 1D multivariate, cacheing point results on other axes
//

double generalIntegratingFunction(double x, void* params) {
	integratingParams* p = (integratingParams*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) {
		if(p->axis < it->second.size()) {
			double r = (it->second)[p->axis];
			assert(r==r);
			return r;
		}
		return 0;
	}
	
	mvec v = p->f(x,p->fparams);
	if(v.size() > p->n_dim) p->n_dim = v.size();
	p->m[x] = v;
	if(p->axis < v.size()) return v[p->axis];
	return 0;
}

/*
 double generalIntegratingFunctionNoCache(double x, void* params) {
 // non-cacheing version always re-evaluates function components
 integratingParams* p = (integratingParams*)params;
 mvec v = p->f(x,p->fparams);
 p->n_dim = v.size();
 assert(p->axis < p->n_dim);
 return v[p->axis];
 }
 */

mvec Integrator::integrate(mvec (*f)(double,void*), double a, double b, void* params) {
	
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

mvec Integrator::_integrate_v(integratingParams& p, double (*integf)(double,void*), double a, double b) {
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
	printf("------------------------------\nGSL Integrator error codes:\n");
	printf("%i: the maximum number of subdivisions was exceeded.\n",GSL_EMAXITER);
    printf("%i: cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.\n",GSL_EROUND);
    printf("%i: a non-integrable singularity or other bad integrand behavior was found in the integration interval.\n",GSL_ESING);
	printf("%i: the integral is divergent, or too slowly convergent to be integrated numerically.\n",GSL_EDIVERGE);
	printf("\n");
}





//
// 2D integrator
//


void Integrator2D::setup_singularities(double a, double b) {
	singularities.clear();
	for(std::vector<vec2>::iterator it = xysingularities.begin(); it != xysingularities.end(); it++)
		singularities.insert((*it)[0]);
	Integrator::setup_singularities(a,b);
}

/// Contains arguments for general 2D integrating function
class integratingParams_2D {
public:
	double (*f2)(vec2,void*);	//< pointer to the 2D function being integrated
	void* fparams;				//< parameters for integrating function
	Integrator* yIntegrator;	//< Integrator for y-direction integrals
	double x;					//< x value being integrated
	double y0;					//< y lower bound of integration
	double y1;					//< y upper bound of integration
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

double Integrator2D::integrate2D(double (*f)(vec2,void*), vec2 ll, vec2 ur, void* params) {
	
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
	double x;					//< x value being integrated
	double y0;					//< y lower bound of integration
	double y1;					//< y upper bound of integration
};

// slice of a 2D function at constant x
mvec xslice_v(double y, void* params) {
	integratingParams_V2D* p = (integratingParams_V2D*)params;
	return p->f2(vec2(p->x,y),p->fparams);
}


double xslice_v_integral(double x, void* params) {
	
	integratingParams_V2D* p = (integratingParams_V2D*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) {
		if(p->axis < it->second.size()) return (it->second)[p->axis];
		return 0;
	}
	
	integratingParams_V2D p2;
	p2.f2 = p->f2;
	p2.fparams = p->fparams;
	p2.x = x;
	p2.yIntegrator = p->yIntegrator;
	p2.y0 = p->y0;
	p2.y1 = p->y1;
	mvec v = p->yIntegrator->integrate(&xslice_v, p->y0, p->y1, &p2);
	
	if(v.size() > p->n_dim) p->n_dim = v.size();
	p->m[x] = v;
	if(p->axis < v.size()) return v[p->axis];
	return 0;
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
	return _integrate_v(p, &xslice_v_integral, ll[0], ur[0]);
}

// polar

/// Contains arguments for polar 2D integrating function
class integratingParams_V2D_P: public integratingParams {
public:
	mvec (*f2)(vec2,void *);	//< pointer to the 2D function being integrated
	Integrator* yIntegrator;	//< Integrator for y-direction integrals
	double r;					//< r distance
	vec2 c;						//< integration center point
	vec2 ll;					//< integration rectangle lower corner
	vec2 ur;					//< integration rectangle upper corner
};

// arc of a 2D function
mvec polar_slice_v(double y, void* params) {
	integratingParams_V2D_P* p = (integratingParams_V2D_P*)params;
	return p->f2( p->c + vec2(cos(y),sin(y))*p->r, p->fparams);
}

angular_interval clip_interval(double l, double r) {
	if(fabs(r) <= fabs(l)) {
		if(l>0) return angular_interval(0,0);
		return angular_interval(0,2*M_PI);
	}
	double a = atan2(sqrt(r*r-l*l),l);
	return angular_interval(-a,a);
}

// determine angular clipping
std::vector<angular_interval> rectangle_clip(vec2 c, vec2 ll, vec2 ur, double r) {
	
 	Angular_Interval_Set AIS;
	AIS.add_interval(0,2*M_PI);
	
	for(unsigned int d=0; d<2; d++) {
		angular_interval a1 = clip_interval(c[d]-ll[d],r);
		a1.add((1+0.5*d)*M_PI);
		angular_interval a2 = clip_interval(ur[d]-c[d],r);
		a2.add(0.5*d*M_PI);
		
		AIS.subtract_interval(a1);
		AIS.subtract_interval(a2);
	}
	
	return  AIS.get_intervals();
}

double polar_slice_v_integral(double x, void* params) {
	
	if(x<=0) return 0;
	assert(params);
	
	integratingParams_V2D_P* p = (integratingParams_V2D_P*)params;
	std::map<double,mvec>::iterator it = p->m.find(x);
	if(it != p->m.end()) {
		if(p->axis < it->second.size()) return (it->second)[p->axis];
		return 0;
	}
	
	//std::cout << "Polar slice at r = " << x << std::endl;
	
	integratingParams_V2D_P p2;
	p2.f2 = p->f2;
	p2.fparams = p->fparams;
	p2.r = x;
	p2.yIntegrator = p->yIntegrator;
	p2.c = p->c;
	
	std::vector<angular_interval> ii = rectangle_clip(p->c, p->ll, p->ur, x);
	if(!ii.size()) return 0;
	
	mvec v;
	for(std::vector<angular_interval>::iterator it = ii.begin(); it != ii.end(); it++) {
		mvec vi = p->yIntegrator->integrate(&polar_slice_v, it->th0, it->th1, &p2);
		if(it==ii.begin()) v = vi;
		else v += vi;
	}
	v *= x;
	
	if(v.size() > p->n_dim) p->n_dim = v.size();
	p->m[x] = v;
	if(p->axis < v.size()) return v[p->axis];
	return 0;
}

double r_max(vec2 ll, vec2 ur, vec2 c) {
	double r[4];
	r[0] = (c-ll).mag2();
	r[1] = (c-vec2(ll[0],ur[1])).mag2();
	r[2] = (c-vec2(ur[0],ll[1])).mag2();
	r[3] = (c-ur).mag2();
	return sqrt(*std::max_element(r,r+4));
}

mvec Integrator2D::polarIntegrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, vec2 c, void* params, double r1, double r0) {
	
	if(r1==-666) r1 = r_max(ll, ur, c); // auto-range
	
	integratingParams_V2D_P p;
	p.yIntegrator = &yIntegrator;
	p.f2 = f;
	p.fparams = params;
	p.ll = ll;
	p.ur = ur;
	p.c = c;
	
	return _integrate_v(p, &polar_slice_v_integral, r0, r1);
}




//
//
//

struct ND_Integrating_Params {
	double (*f1)(mvec, void*);
	mvec (*fN)(mvec, void*);
	void* fparams;
};

// function wrapper for cubature
int cubature_f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	ND_Integrating_Params* p = (ND_Integrating_Params*)fdata;
	mvec vx(x, x+ndim);
	if(p->fN) {
		mvec vy = (*p->fN)(vx, p->fparams);
		for(unsigned int i=0; i<fdim; i++) fval[i] = vy[i];
	} else if(p->f1) {
		assert(fdim == 1);
		fval[0] = (*p->f1)(vx, p->fparams);
	} else return -1;
    return 0;
}

IntegratorND::IntegratorND(): rel_err(1e-4), abs_err(1e-5) {}

double IntegratorND::integrate(double (*f)(mvec,void*), mvec ll, mvec ur, void* params) const {
	ND_Integrating_Params p;
	p.fN = NULL;
	p.f1 = f;
	p.fparams = params;
	
	assert(ll.size() == ur.size());
	
	double val;
	double err;
	int i = hcubature(1, &cubature_f, &p,
					  ll.size(), &ll[0], &ur[0],
					  0, abs_err, rel_err, ERROR_L2,
					  &val, &err);
	if(i) printf("(*Integration Warning: cubature error %i*)\n",i);
	return val;
}

mvec IntegratorND::integrate(mvec (*f)(mvec,void*), unsigned int fdim, mvec ll, mvec ur, void* params) const {
	ND_Integrating_Params p;
	p.fN = f;
	p.f1 = NULL;
	p.fparams = params;
	
	assert(ll.size() == ur.size());
	
	mvec val(fdim);
	mvec err(fdim);
	int i = hcubature(fdim, &cubature_f, &p,
					  ll.size(), &ll[0], &ur[0],
					  0, abs_err, rel_err, ERROR_L2,
					  &val[0], &err[0]);
	if(i) printf("(*Integration Warning: cubature error %i*)\n",i);
	return val;
}


