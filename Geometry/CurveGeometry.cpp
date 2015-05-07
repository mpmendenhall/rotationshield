#include "CurveGeometry.hh"

CurveGeometry::CurveGeometry(): dflt_integrator_ndivs(5) {
    myIntegrator.setMethod(INTEG_GSL_QAG);
}

double integration_dl(double x, void* params) {
    CurveGeometry* S = (CurveGeometry*)params;
    return S->dl(x);
}

double CurveGeometry::length(double x0, double x1) const {
    return myIntegrator.integrate(&integration_dl, x0, x1, (void*)this);
}

mvec CurveGeometry::subdividedIntegral(mvec (*f)(double, void*), void* fparams, double x0, double x1, unsigned int ndx) const {
    ndx = ndx?ndx:dflt_integrator_ndivs;
    x1 = (x1-x0)/ndx;
  
    mvec m;
    for(unsigned int nx=0; nx<ndx; nx++) {
        double xx = x0 + x1*nx;
        mvec mi = myIntegrator.integrate(f, xx, xx+x1, fparams);
        if(!nx) m = mi;
        else m += mi;
    }
    return m;
}

//////////////////////////////////////////////

vec3 EllipseCurve::operator()(double x) const {
    return x0 + a1*cos(2*M_PI*x) + a2*sin(2*M_PI*x);
}

CircleCurve::CircleCurve(vec3 xx0, vec3 dx): EllipseCurve(xx0, vec3(), vec3()), r(dx.mag()) {
    /// TODO calculate axes
}
