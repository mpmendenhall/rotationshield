#include "SurfaceGeometry.hh"
#include "Integrator.hh"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>

SurfaceGeometry::SurfaceGeometry(): dflt_integrator_ndivs_x(8), dflt_integrator_ndivs_y(8),
polar_integral_center(NULL), polar_r0(0) {
    myIntegrator.setMethod(INTEG_GSL_QAG);
    myIntegrator2D.setMethod(INTEG_GSL_QAG);
    myMinimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, 2);
    assert(myMinimizer);
}

SurfaceGeometry::~SurfaceGeometry() {
    gsl_multimin_fdfminimizer_free(myMinimizer);
}

vec3 SurfaceGeometry::snorm(const vec2& p, bool normalized) const {
    return cross(deriv(p,0,normalized), deriv(p,1,normalized));
}

vec2 SurfaceGeometry::d_pathlength(vec2 l) const {
    return vec2(deriv(l,0).mag(), deriv(l,1).mag());
}

Matrix<3,3,double> SurfaceGeometry::rotToLocal(const vec2& x) const {
    vec3 v0 = deriv(x,0,true);
    vec3 v1 = deriv(x,1,true);
    vec3 v2 = cross(v0,v1);
    Matrix<3,3,double> M;
    for(unsigned int i=0; i<3; i++) {
        M(0,i) = v0[i];
        M(1,i) = v1[i];
        M(2,i) = v2[i];
    }
    return M;
}

double SurfaceGeometry::dA(const vec2& l) const {
    return snorm(l,false).mag();
}

double integration_dA(mvec l, void* params) {
    SurfaceGeometry* S = (SurfaceGeometry*)params;
    return S->dA(vec2(l[0],l[1]));
}

double SurfaceGeometry::area(const vec2& ll, const vec2& ur) const {
    return myIntegratorND.integrate(&integration_dA, ll, ur, (void*)this);
}

void SurfaceGeometry::cache_sincos(double theta, double& s, double& c) const {
    static double c_th = theta;
    static double c_s = sin(theta);
    static double c_c = cos(theta);
    if(c_th != theta) { c_th=theta; c_c=cos(c_th); c_s=sin(c_th); }
    s = c_s;
    c = c_c;
}

void SurfaceGeometry::proximity(vec3 p, vec2 ll, vec2 ur, double& mn, double& mx) const {
    const int nx = 4;
    const int ny = 4;
    double r[nx*ny];
    
    unsigned int i=0;
    for(int x = 0; x<nx; x++) {
        for(int y = 0; y<ny; y++) {
            double l1 = double(x)/(nx-1);
            double l2 = double(y)/(ny-1);
            vec2 v( ll[0]*(1-l1)+ur[0]*l1, ll[1]*(1-l2)+ur[1]*l2);
            r[i++] = ( p - (*this)(v) ).mag2();
        }
    }
    
    mn = sqrt(*std::max_element(r,r+nx*ny));
    mx = sqrt(*std::min_element(r,r+nx*ny));
}

struct f2_to_fN_params {
    mvec (*f2)(vec2, void*);
    void* fparams;
};

mvec f2_to_fN(mvec v, void* pp) {
    f2_to_fN_params* p = (f2_to_fN_params*)pp;
    return (*p->f2)(vec2(v[0],v[1]), p->fparams);
}

mvec SurfaceGeometry::subdividedIntegral(mvec (*f)(vec2, void*), unsigned int fdim, void* fparams, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
    
    ndx = ndx?ndx:dflt_integrator_ndivs_x;
    ndy = ndy?ndx:dflt_integrator_ndivs_y;
    for(unsigned int i=0; i<2; i++)
        ur[i] = (ur[i]-ll[i])/(i?ndy:ndx);
    
    mvec m;
    for(unsigned int nx=0; nx<ndx; nx++) {
        for(unsigned int ny=0; ny<ndy; ny++) {
            vec2 lll = ll + ur*vec2(nx,ny);
            mvec mi;
            if(polar_integral_center) {
                mi = myIntegratorND.integratePolar(f, fdim, *polar_integral_center, lll, lll+ur, fparams, -666, polar_r0);
            } else if(myIntegrator2D.getMethod() == INTEG_GSL_CQUAD) {
                mi = myIntegrator2D.integrate2D(f, lll, lll+ur, fparams);
            } else {
                f2_to_fN_params p;
                p.f2 = f;
                p.fparams = fparams;
                mi = myIntegratorND.integrate(&f2_to_fN, fdim, mvec(lll), mvec(lll+ur), &p);
            }
            if(!nx && !ny) m = mi;
            else m += mi;
        }
    }
    return m;
}

struct closestPoint_params {
    const SurfaceGeometry* S;
    vec3 x;
};

double surface_distance_f(const gsl_vector* l, void * params) {
    closestPoint_params* p = (closestPoint_params*)params;
    return (p->x - (*p->S)(vec2(gsl_vector_get(l,0),gsl_vector_get(l,1)))).mag2();
}

void surface_distance_df(const gsl_vector* x, void * params, gsl_vector * g) {
    closestPoint_params* p = (closestPoint_params*)params;
    vec2 l(gsl_vector_get(x,0),gsl_vector_get(x,1));
    vec3 r = p->x - (*p->S)(l);
    for(unsigned int i=0; i<2; i++)
        gsl_vector_set(g,i,-2*r.dot(p->S->deriv(l,i)));
}

void surface_distance_fdf(const gsl_vector* l, void * params, double * f, gsl_vector * g) {
    *f = surface_distance_f(l, params);
    surface_distance_df(l, params, g);
}

vec2 SurfaceGeometry::closestPoint(vec3 x, double& d2) const {
    
    // initial guess
    unsigned int nx = 5;
    unsigned int ny = 5;
    gsl_vector* l = gsl_vector_alloc(2);
    d2 = DBL_MAX;
    for(unsigned int ix = 0; ix<nx; ix++) {
        for(unsigned int iy = 0; iy<ny; iy++) {
            vec2 l0((ix+0.5)/double(nx),(iy+0.5)/double(ny));
            double d2i = (x-(*this)(l0)).mag2();
            if(d2i < d2) {
                d2 = d2i;
                gsl_vector_set(l, 0, l0[0]);
                gsl_vector_set(l, 1, l0[1]);
            }
        }
    }
    
    // minimizer setup
    closestPoint_params p;
    p.S = this;
    p.x = x;
    gsl_multimin_function_fdf my_func;
    my_func.n = 2;
    my_func.f = &surface_distance_f;
    my_func.df = &surface_distance_df;
    my_func.fdf = &surface_distance_fdf;
    my_func.params = &p;
    gsl_multimin_fdfminimizer_set(myMinimizer, &my_func, l, 0.01, 0.1);
    
    // iterate to minimum
    int status;
    do {
        status = gsl_multimin_fdfminimizer_iterate(myMinimizer);
        if(status) break;
        status = gsl_multimin_test_gradient(myMinimizer->gradient, 1e-6);
    } while(status == GSL_CONTINUE);
    
    // cleanup and return
    vec2 l0(gsl_vector_get(gsl_multimin_fdfminimizer_x(myMinimizer), 0), gsl_vector_get(gsl_multimin_fdfminimizer_x(myMinimizer), 1));
    d2 =  gsl_multimin_fdfminimizer_minimum(myMinimizer);
    gsl_vector_free(l);
    return l0;
}


//--------------------------------------

vec2 CylSurfaceGeometry::cache_profile(double l) const {
    assert(zr_profile);
    static double c_l = l;
    static vec2 c_v = (*zr_profile)(l);
    if(c_l != l) { c_l = l; c_v = (*zr_profile)(l); }
    return c_v;
}

vec3 CylSurfaceGeometry::operator()(const vec2& p) const {
    vec2 zr = cache_profile(p[0]);
    double phi = 2*M_PI*p[1];
    double s,c;
    cache_sincos(phi,s,c);
    return vec3(zr[1]*c, zr[1]*s, zr[0]);
}

vec3 CylSurfaceGeometry::deriv(const vec2& p, unsigned int i,  bool normalized) const {
    
    double phi = 2*M_PI*p[1];
    double s,c;
    cache_sincos(phi,s,c);
    assert(zr_profile);
    
    if(i==0) {
        vec2 dzr = zr_profile->deriv(p[0],normalized);
        assert(dzr.mag2());
        return vec3(dzr[1]*c, dzr[1]*s, dzr[0]);
    } else if(i==1) {
        vec2 zr = cache_profile(p[0]);
        if(normalized) return vec3(-s*sign(zr[1]), c*sign(zr[1]), 0);
        return vec3(-zr[1]*s*2*M_PI, zr[1]*c*2*M_PI, 0);
    }
    
    assert(false);
    return vec3(0,0,0);
}

vec3 CylSurfaceGeometry::snorm(const vec2& p, bool normalized) const {
    if(!normalized) return SurfaceGeometry::snorm(p, false);
    
    double phi = 2*M_PI*p[1];
    double s,c;
    cache_sincos(phi,s,c);
    vec2 dzr = zr_profile->deriv(p[0],true);
    return vec3(-dzr[0]*c, -dzr[0]*s, dzr[1]);
}

double CylSurfaceGeometry::dA(const vec2& l) const {
    vec2 zr = cache_profile(l[0]);
    vec2 dzr = zr_profile->deriv(l[0]);
    return 2*M_PI*zr[1]*dzr.mag();
}

double integration_cyl_dA(double x, void* params) {
    CylSurfaceGeometry* S = (CylSurfaceGeometry*)params;
    return S->dA(vec2(x,0));
}

double CylSurfaceGeometry::area(const vec2& ll, const vec2& ur) const {
    return myIntegrator.integrate(&integration_cyl_dA, ll[0], ur[0], (void*)this)*(ur[1]-ll[1]);
}
