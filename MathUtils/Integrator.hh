/// \file integrator.hh
/// \brief 3-vector integration
///
/// Extends the GSL numerical integration routines to work
/// with vector-valued functions, by integrating over each
/// component separately. The 3-vector values evaluated at
/// each point are cached so that they can be used by all
/// three component integrals.

#ifndef INTEGRATOR_HH
/// Make sure the file is only loaded once
#define INTEGRATOR_HH 1

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "VarVec.hh"
#include "Vec.hh"
#include <map>


class Integrator;

/// Contains the arguments for generalIntegratingFunction()
struct integratingParams {
	void* fparams; //< additional arguments for the function being integrated
	vec3 (*f)(mdouble,void *); //< pointer to the function being integrated
	int axis; //< which of the three vector components is currently being integrated
	std::map<double,vec3> m; //< cache for evaluated function points
	bool verbose; //< whether to display pts during integration
};

/// Used as a wrapper for 3-vector functions, presenting the correct form needed by the GSL integration routines
double generalIntegratingFunction(double x, void* params);

///	Integration of vector-valued functions

///	The Integrator class uses gsl_integration_qag()
///	running in GSL_INTEG_GAUSS15 mode to evaluate
///	vector-valued functions, caching the evaluated
///	points in a SequenceCache so they can be used
///	for all three component integrals.
class Integrator {
public:
	/// Constructor
	Integrator() { gslIntegrationWS = gsl_integration_workspace_alloc(512); gsl_set_error_handler_off(); }
	/// Destructor
	~Integrator() { gsl_integration_workspace_free(gslIntegrationWS); }
	/// Integrates a vector-valued function \f$ \int_a^b \vec f(x)dx\f$ using GSL numerical integration routines for each component
	/** \param f vector-valued function of a real variable to be integrated
	 \param a lower bound of integration
	 \param b upper bound of integration
	 \param params additional parameters for the integrated function */
	vec3 integrate(vec3 (*f)(mdouble,void*),mdouble a, mdouble b, void* params = 0x0) {
		vec3 v; integratingParams p;
		p.fparams = params;
		p.f = f;
		p.m = std::map<double,vec3>();
		gsl_function F;
		F.function = &generalIntegratingFunction;
		F.params = &p;
		double r,e;
		p.verbose = false;
		size_t neval;
		
		for(int i=0; i<3; i++) {
			p.axis = i;
			int er = gsl_integration_qng(&F, a, b, 1e-7, 1e-7, &r, &e, &neval);
			if(er)
				printf("(*INTEGRATION WARNING*)\n");
			v[i] = r;
		}
		return v;
	}
	
	gsl_integration_workspace* gslIntegrationWS; //< needed by GSL integration routines called in integrate()
};


#endif
