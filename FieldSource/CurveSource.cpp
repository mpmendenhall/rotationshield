#include "CurveSource.hh"
#include <cassert>

CurveSource::CurveSource(CurveGeometry* CG, double jj, const string& nm): FieldSource(nm), myCurve(CG), vis_n(36), j0(jj) {
    if(myCurve) mySymmetry = myCurve->getSymmetry();
}

vec3 CurveSource::fieldAt_contrib_from(const vec3& v, double x) const {
    // Biot-Savart law
    vec3 x0 = (*myCurve)(x);
    vec3 r = v - x0;
    double mr = r.mag();
    if(!mr) return vec3(0,0,0);
    vec3 dI = myCurve->deriv(x)*j0;
    vec3 B = cross(dI,r)/(4.*M_PI*mr*mr*mr);
    assert(B[0]==B[0] && B[1]==B[1] && B[2]==B[2]);
    return B;
}

struct CurveSourceIntegParams {
    const CurveSource* S;
    vec3 v;
    const Matrix<2,3,double>* M2;
    const Matrix<3,3,double>* M3;
};

mvec CSdl(double x, void* params) {
    CurveSourceIntegParams& p = *(CurveSourceIntegParams*)params;
    vec3 B = p.S->fieldAt_contrib_from(p.v,x);
    if(p.M2) return mvec( (*p.M2)*B );
    if(p.M3) return mvec( (*p.M3)*B );
    return mvec(B);
}

vec3 CurveSource::fieldAt(const vec3& v, double x0, double x1, unsigned int ndx) const {
    CurveSourceIntegParams p;
    p.S = this;
    p.v = v;
    p.M2 = NULL;
    p.M3 = NULL;
    assert(myCurve);
    
    mvec B = myCurve->subdividedIntegral(&CSdl, &p, x0, x1, ndx);
    assert(B.size()==3);
    return vec3(B[0],B[1],B[2]);
}

vec2 CurveSource::fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M, double x0, double x1, unsigned int ndx) const {
    CurveSourceIntegParams p;
    p.S = this;
    p.v = v;
    p.M2 = &M;
    p.M3 = NULL;
    assert(myCurve);
    
    mvec MB = myCurve->subdividedIntegral(&CSdl, &p, x0, x1, ndx);
    assert(MB.size()==2);
    return vec2(MB[0],MB[1]);
}

vec3 CurveSource::fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M, double x0, double x1, unsigned int ndx) const {
    CurveSourceIntegParams p;
    p.S = this;
    p.v = v;
    p.M2 = NULL;
    p.M3 = &M;
    assert(myCurve);
    
    mvec MB = myCurve->subdividedIntegral(&CSdl, &p, x0, x1, ndx);
    assert(MB.size()==3);
    return vec3(MB[0],MB[1],MB[2]);
}

void CurveSource::_visualize() const {
    assert(myCurve);
    vsr::startLines();
    for(unsigned int i=0; i<vis_n; i++) {
        float p = i/double(vis_n-1);
        float q = 1.0-p;
        if(p<0.5)
            vsr::setColor(2*p,0,1.0,1.0);
        else
            vsr::setColor(1.0,0,2*q,1.0);
        vsr::vertex((*myCurve)(j0 >= 0? p : q));
    }
    vsr::endLines();
}
