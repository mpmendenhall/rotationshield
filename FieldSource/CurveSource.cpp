#include "CurveSource.hh"

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
    
    mvec MB = myCurve->subdividedIntegral(&CSdl, &p, x0, x1, ndx);
    assert(MB.size()==3);
    return vec3(MB[0],MB[1],MB[2]);
}
