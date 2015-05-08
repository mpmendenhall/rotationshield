#ifndef CURVEGEOMETRY_HH
/// Make sure this header is only loaded once
#define CURVEGEOMETRY_HH

#include "Typedefs.hh"
#include "DVFunc.hh"
#include "Integrator.hh"
#include "SymmetrySpec.hh"

/// Parametric curve in 3D space
class CurveGeometry: public DVFunc1<3,double> {
public:
    /// Constructor
    CurveGeometry();
    /// destructor
    virtual ~CurveGeometry() { }
        
    /// evaluate function
    virtual vec3 operator()(double x) const = 0;
    /// evaluate differential pathlength at position
    virtual double dl(double x) const { return deriv(x).mag(); }
    
    /// length corresponding to x range
    virtual double length(double x0 = 0, double x1 = 1) const;
        
    /// whether curve is closed
    virtual bool isClosed() const { return false; }
    
    /// return information on surface symmetries
    const SymmetrySpec& getSymmetry() const { return mySymmetry; }
        
    unsigned int dflt_integrator_ndivs; ///< default number of sections to partition x integral in
        
    /// convenience mechanism for line integrals
    mvec subdividedIntegral(mvec (*f)(double, void*), void* fparams, double x0, double x1, unsigned int ndx = 0) const;
        
    mutable Integrator myIntegrator;    ///< integrator for internal calculations
        
protected:
    SymmetrySpec mySymmetry;            ///< symmetries of curve
};


/// Ellipse curve
class EllipseCurve: public CurveGeometry {
public:
    /// Constructor
    EllipseCurve(vec3 xx0, vec3 aa1, vec3 aa2): x0(xx0), a1(aa1), a2(aa2) { }
    
    /// evaluate function
    virtual vec3 operator()(double x) const;
    /// derivative
    virtual vec3 deriv(double x, bool normalized = false) const;
    /// whether curve is closed
    virtual bool isClosed() const { return true; }
    
protected:
    vec3 x0;            ///< center
    vec3 a1;            ///< first axis (position at x=0)
    vec3 a2;            ///< second axis (position at x=1/4)
};

/// Circle, as special case of ellipse
class CircleCurve: public EllipseCurve {
public:
    /// Constructor
    CircleCurve(vec3 xx0, vec3 dx);
    /// evaluate differential pathlength at position
    virtual double dl(double) const { return 2*M_PI*r; }
    /// length corresponding to x range
    virtual double length(double x0 = 0, double x1 = 1) const { return (x1-x0)*2*M_PI*r; }
    
protected:
    double r;           ///< radius
};

#endif
