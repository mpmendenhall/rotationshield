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
#include "Typedefs.hh"
#include <map>

///	Integration of vector-valued functions

/// Contains the arguments for generalIntegratingFunction()
class integratingParams {
public:
	/// destructor
	virtual ~integratingParams() {}
	
	void* fparams; 					//< additional arguments for the function being integrated
	mvec (*f)(mdouble,void *);		//< pointer to the function being integrated
	unsigned int axis;				//< which of the three vector components is currently being integrated
	unsigned int n_dim;				//< dimension of return vector
	std::map<double,mvec> m;		//< cache for evaluated function points
	
	static bool verbose;			//< whether to display pts during integration
};

///	The Integrator class uses gsl_integration_qag()
///	running in GSL_INTEG_GAUSS15 mode to evaluate
///	vector-valued functions, caching the evaluated
///	points so they can be used for all three component integrals.
class Integrator {
public:
	/// Constructor
	Integrator(): gslIntegrationWS(gsl_integration_workspace_alloc(512)) { gsl_set_error_handler_off(); }
	/// Destructor
	virtual ~Integrator() { gsl_integration_workspace_free(gslIntegrationWS); }
	/// Integrates a vector-valued function \f$ \int_a^b \vec f(x)dx\f$ using GSL numerical integration routines for each component
	/** \param f vector-valued function of a real variable to be integrated
	 \param a lower bound of integration
	 \param b upper bound of integration
	 \param params additional parameters for the integrated function */
	mvec integrate(mvec (*f)(mdouble,void*), mdouble a, mdouble b, void* params = 0x0);

protected:
	/// internal calls to GSL integration
	mvec _integrate(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b);
	
	gsl_integration_workspace* gslIntegrationWS; //< needed by GSL integration routines called in integrate()
};

/// 2-dimensional vector integrator
class Integrator2D: public Integrator {
public:
	/// Constructor
	Integrator2D(): Integrator() {}

	/// Integrates a vector-valued function \f$ \int_{x_0}^{x_1} \int_{y_0}^{y_1} \vec f(x,y) dy dx\f$
	/** \param f vector-valued function of a real variable to be integrated
	 \param x0 lower x bound of integration
	 \param x1 upper x bound of integration
	 \param y0 lower y bound of integration
	 \param y1 upper y bound of integration
	 \param params additional parameters for the integrated function */
	mvec integrate(mvec (*f)(mdouble,mdouble,void*), mdouble x0, mdouble x1, mdouble y0, mdouble y1, void* params = 0x0);
};



#endif
