#include "CurveGeometry.hh"

CurveGeometry::CurveGeometry(): dflt_integrator_ndivs(1) {
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

vec3 EllipseCurve::deriv(double x, bool normalized) const {
    vec3 dx = (a1*(-sin(2*M_PI*x)) + a2*cos(2*M_PI*x))*2*M_PI;
    if(normalized) return dx.normalized();
    return dx;
}

CircleCurve::CircleCurve(vec3 xx0, vec3 dx): EllipseCurve(xx0, vec3(), vec3()), r(dx.mag()) {
    double m0 = fabs(dx[0]);
    double m1 = fabs(dx[1]);
    double m2 = fabs(dx[2]);
    
    // TODO check consistency of cross-products signs
    if(m0 > m1) {
        if(m1 > m2) a1 = vec3(-dx[1],dx[0],0);
        else a1 = vec3(-dx[2],0,dx[0]);
    } else {
        if(m0 > m2) a1 = vec3(-dx[1],dx[0],0);
        else a1 = vec3(0,-dx[2],dx[1]);
    }
    
    dx /= r;
    a1 = a1.normalized()*r;
    a2 = cross(dx, a1);
    
    if(!dx[0] && !dx[1]) mySymmetry.rotation = true;
}
