/// \file integrator.hh
/// \brief 3-vector integration
///
/// Extends the GSL numerical integration routines to work
/// with vector-valued functions, by integrating over each
/// component separately. The 3-vector values evaluated at
/// each point are cached so that they can be used by all
/// three component integrals.

/* 
 * Integrator.hh, part of the RotationShield program
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

#ifndef INTEGRATOR_HH
/// Make sure this header is only loaded once
#define INTEGRATOR_HH

#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "Typedefs.hh"
#include <map>
#include <string>
#include <set>

///	Integration of vector-valued functions

/// internal integration method specifier
enum Integration_Method {
	INTEG_GSL_QNG  = 0,	///< non-adaptive
	INTEG_GSL_QAG  = 1,	///< adaptive
	INTEG_GSL_QAGS = 2,	///< adaptive with singular points
	INTEG_GSL_QAGP = 3,	///< adaptive with specified singular points
	INTEG_GSL_CQUAD= 4	///< adaptive CQUAD
};

const unsigned int INTEG_WS_SIZE = 512;	///< size of integrating workspace (max. number of evaluation intervals)

/// Contains the arguments for generalIntegratingFunction()
class integratingParams {
public:
	/// constructor
	integratingParams(): axis(0), n_dim(0) {}
	
	/// destructor
	virtual ~integratingParams() {}
	
	void* fparams; 					///< additional arguments for the function being integrated
	mvec (*f)(double,void *);		///< pointer to the function being integrated
	unsigned int axis;				///< which of the three vector components is currently being integrated
	unsigned int n_dim;				///< dimension of return vector
	std::map<double,mvec> m;		///< cache for evaluated function points
	
	std::string nm;					///< name describing integration, for error display
	static bool verbose;			///< whether to display pts during integration
};

///	The Integrator class uses gsl_integration_qag()
///	running in GSL_INTEG_GAUSS15 mode to evaluate
///	vector-valued functions, caching the evaluated
///	points so they can be used for all three component integrals.
class Integrator {
public:
	/// Constructor
	Integrator();
	/// Destructor
	virtual ~Integrator();
	
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
	mvec integrate(mvec (*f)(double,void*), double a, double b, void* params = NULL);
	
	double rel_err;					///< relative error target, OR
	double abs_err;					///< absolute error target
	unsigned int err_count;			///< note of integration errors
	/// return and reset error counter
	unsigned int reset_errcount() { unsigned int e = err_count; err_count = 0; return e; }
	
	std::set<double> singularities;	///< known singular points
	
	/// print meanings of GSL error codes
	static void printErrorCodes();
	
	/// get integration method
	Integration_Method getMethod() const { return myMethod; }
	/// set integration method
	virtual void setMethod(Integration_Method m) { myMethod = m; }
	
protected:
	/// internal multi-variate integration
	mvec _integrate_v(integratingParams& p, double (*integf)(double,void*), double a, double b);
	
	/// internal calls to GSL integration
	double _integrate(gsl_function* F, double a, double b);
	
	/// set list of interesting singularities
	virtual void setup_singularities(double a, double b);
	std::vector<double> _singularities;	///< singularities in integrating range
	
	Integration_Method myMethod;					///< selection of internal integration method
	gsl_integration_workspace* gslIntegrationWS; 	///< needed by GSL integration routines called in integrate()
	gsl_integration_cquad_workspace * gsl_cqd_ws;	///< needed by
};

/// 2-dimensional vector integrator; main use is polar integration around singularity (otherwise, use IntegratorND below)
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
	double integrate2D(double (*f)(vec2,void*), vec2 ll, vec2 ur, void* params = NULL);
	
	/// Integrates a vector-valued function \f$ \int_{x_0}^{x_1} \int_{y_0}^{y_1} \vec f(x,y) dy dx\f$
	/** \param f vector-valued function of a real variable to be integrated
	 \param x0 lower x bound of integration
	 \param x1 upper x bound of integration
	 \param y0 lower y bound of integration
	 \param y1 upper y bound of integration
	 \param params additional parameters for the integrated function */
	mvec integrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, void* params = NULL);
	
	/// performs polar integral around point, clipped to specfied rectangle
	mvec polarIntegrate2D(mvec (*f)(vec2,void*), vec2 ll, vec2 ur, vec2 c, void* params = NULL, double r1 = -666, double r0 = 0);
	
	std::vector<vec2> xysingularities;	///< known singularities
	
	/// set integration method
	virtual void setMethod(Integration_Method m) { Integrator::setMethod(m); yIntegrator.setMethod(m); }
	
	/// return and reset error counter in y direction
	unsigned int reset_y_errcount() { return yIntegrator.reset_errcount(); }
	
protected:
	
	/// set list of interesting singularities
	virtual void setup_singularities(double a, double b);
	
	Integrator yIntegrator; ///< integrator for second dimension
};

/// Multi-dimensional integrator, wrapper for "cubature" code
class IntegratorND {
public:
	/// Constructor
	IntegratorND();
	/// integrate scalar-valued function
	double integrate(double (*f)(mvec,void*), mvec ll, mvec ur, void* params = NULL) const;
	/// integrate vector-valued function
	mvec integrate(mvec (*f)(mvec,void*), unsigned int fdim, mvec ll, mvec ur, void* params = NULL) const;
	/// perform 2D integral in polar form, with internal sum of opposite angles
	mvec integratePolar(mvec (*f)(vec2, void*), unsigned int fdim, vec2 x0, vec2 ll, vec2 ur, void* params = NULL, double r1 = -666, double r0 = 0) const;
	
	
	double rel_err;		///< relative error target, OR
	double abs_err;		///< absolute error target
};


#endif
