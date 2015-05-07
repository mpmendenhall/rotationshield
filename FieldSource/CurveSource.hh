/// \file CurveSource.hh \brief Curved current line sources

#ifndef CURVESOURCE_HH
/// Make sure this header is only loaded once
#define CURVESOURCE_HH

#include "FieldSource.hh"
#include "CurveGeometry.hh"

class CurveSource: public FieldSource {
public:
    /// Constructor
    CurveSource(CurveGeometry* CG = NULL, double jj = 0, const string& nm = "CurveSource");
    /// Destructor
    virtual ~CurveSource() { if(ownsCurve) delete myCurve; }
    
    /// field contribution f(x)dl; x in [0,1]
    virtual vec3 fieldAt_contrib_from(const vec3& v, double x) const;
    
    /// Magnetic field at a certain point from a sub-domain
    virtual vec3 fieldAt(const vec3& v, double x0, double x1, unsigned int ndx = 0) const;
    /// Magnetic field with transform at a certain point from a sub-domain
    virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M, double x0, double x1, unsigned int ndx = 0) const;
    /// Magnetic field with transform at a certain point from a sub-domain
    virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M, double x0, double x1, unsigned int ndx = 0) const;
    
    /// Magnetic field at a specified point
    virtual vec3 fieldAt(const vec3& v) const { return fieldAt(v,0,1); }
    /// Magnetic field at a point, with interaction matrix
    virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M) const { return fieldAtWithTransform2(v,M,0,1); }
    /// Magnetic field at a point, with interaction matrix
    virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M) const { return fieldAtWithTransform3(v,M,0,1); }
    
    /// Visualize the field source
    virtual void _visualize() const;
    
    CurveGeometry* myCurve;     ///< surface on which source is defined
    bool ownsCurve = false;     ///< whether this object "owns" myCurve for deletion
    unsigned int vis_n;         ///< visualization divisions
    double j0;                  ///< current along curve
};

#endif
