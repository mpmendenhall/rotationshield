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
#include <string>
#include <set>

///	Integration of vector-valued functions

/// internal integration method specifier
enum Integration_Method {
	INTEG_GSL_QNG  = 0,	//< non-adaptive
	INTEG_GSL_QAG  = 1,	//< adaptive
	INTEG_GSL_QAGS = 2,	//< adaptive with singular points
	INTEG_GSL_QAGP = 3,	//< adaptive with specified singular points
	INTEG_GSL_CQUAD= 4	//< adaptive CQUAD
};

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
	
	std::string nm;					//< name describing integration, for error display
	static bool verbose;			//< whether to display pts during integration
};

///	The Integrator class uses gsl_integration_qag()
///	running in GSL_INTEG_GAUSS15 mode to evaluate
///	vector-valued functions, caching the evaluated
///	points so they can be used for all three component integrals.
class Integrator {
public:
	/// Constructor
	Integrator(): rel_err(1e-4), abs_err(1e-5), myMethod(INTEG_GSL_QNG),
		gslIntegrationWS(gsl_integration_workspace_alloc(512)),
		gsl_cqd_ws(gsl_integration_cquad_workspace_alloc(512)) { gsl_set_error_handler_off(); }
	/// Destructor
	virtual ~Integrator() {
		gsl_integration_workspace_free(gslIntegrationWS);
		gsl_integration_cquad_workspace_free(gsl_cqd_ws);
	}
	
	/// Integrates a function \f$ \int_a^b f(x)dx\f$ using GSL numerical integration routines for each component
	/** \param f function to be integrated
	 \param a lower bound of integration
	 \param b upper bound of integration
	 \param params additional parameters for the integrated function */
	double integrate(double (*f)(double, void*), double a, double b, void* params = NULL);
	
	/// Integrates a vector-valued function \f$ \int_a^b \vec f(x)dx\f$ using GSL numerical integration routines for each component
	/** \param f vector-valued function of a real variable to be integrated
	 \param a lower bound of integration
	 \param b upper bound of integration
	 \param params additional parameters for the integrated function */
	mvec integrate(mvec (*f)(mdouble,void*), mdouble a, mdouble b, void* params = NULL);
	
	double rel_err;					//< relative error target, OR
	double abs_err;					//< absolute error target
	
	std::set<double> singularities;	//< known singular points
	
	/// print meanings of GSL error codes
	static void printErrorCodes();
	
	/// get integration method
	Integration_Method getMethod() const { return myMethod; }
	/// set integration method
	virtual void setMethod(Integration_Method m) { myMethod = m; }
	
protected:
	/// internal multi-variate integration
	mvec _integrate_v(integratingParams& p, double (*integf)(double,void*), mdouble a, mdouble b);
	
	/// internal calls to GSL integration
	double _integrate(gsl_function* F, double a, double b);
	
	/// set list of interesting singularities
	virtual void setup_singularities(mdouble a, mdouble b);
	std::vector<double> _singularities;	//< singularities in integrating range
	
	Integration_Method myMethod;					//< selection of internal integration method
	gsl_integration_workspace* gslIntegrationWS; 	//< needed by GSL integration routines called in integrate()
	gsl_integration_cquad_workspace * gsl_cqd_ws;	//< needed by
};

/// 2-dimensional vector integrator
class Integrator2D: public Integrator {
public:
	/// Constructor
	Integrator2D(): Integrator() {}

	/// Integrates a 2-variable function \f$ \int_{x_0}^{x_1} \int_{y_0}^{y_1} f(x,y) dy dx\f$
	/** \param f vector-valued function of a real variable to be integrated
	 \param x0 lower x bound of integration
	 \param x1 upper x bound of integration
	 \param y0 lower y bound of integration
	 \param y1 upper y bound of integration
	 \param params additional parameters for the integrated function */
	mdouble integrate2D(mdouble (*f)(vec2,void*), vec2 ll, vec2 ur, void* params = NULL);
	
	/// Integrates a vector-valued function \f$ \int_{x_0}^{x_1} \int_{y_0}^{y_1} \vec f(x,y) dy dx\f$
	/** \param f vector-valued function of a real variable to be integrated
	 \param x0 lower x bound of integration
	 \param x1 upper x bound of integration
	 \param y0 lower y bound of integration
	 \param y1 upper y bound of integration
	 \param params additional parameters for the integrated function */
	mvec integrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* params = NULL);
	
	std::vector<vec2> xysingularities;	//< known singularities
	
	/// set integration method
	virtual void setMethod(Integration_Method m) { Integrator::setMethod(m); yIntegrator.setMethod(m); }
	
protected:
	
	/// set list of interesting singularities
	virtual void setup_singularities(mdouble a, mdouble b);
	
	Integrator yIntegrator; //< integrator for second dimension
};



#endif
