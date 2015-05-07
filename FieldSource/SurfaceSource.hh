/// \file SurfaceSource.hh \brief Base class for extended 2D surface magnetic field sources

#ifndef SURFACESOURCE_HH
/// Make sure this header is only loaded once
#define SURFACESOURCE_HH

#include "FieldSource.hh"
#include "Integrator.hh"
#include "SurfaceGeometry.hh"

class SurfaceSource: public FieldSource {
public:
    /// constructor
    SurfaceSource(SurfaceGeometry* SG = NULL, const string& nm = "SurfaceSource"):
        FieldSource(nm), mySurface(SG), vis_n1(200), vis_n2(200) {}

    /// destructor
    virtual ~SurfaceSource() {}
    
    /// field contribution f(x,y)dA; x,y in [0,1]^2
    virtual vec3 fieldAt_contrib_from(const vec3& v, const vec2& l) const = 0;
    
    /// Magnetic field at a certain point from a sub-domain
    virtual vec3 fieldAt(const vec3& v, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
    /// Magnetic field with transform at a certain point from a sub-domain
    virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
    /// Magnetic field with transform at a certain point from a sub-domain
    virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M, vec2 ll, vec2 ur, unsigned int ndx = 0, unsigned int ndy = 0) const;
    
    /// Magnetic field at a specified point
    virtual vec3 fieldAt(const vec3& v) const { return fieldAt(v,vec2(0,0),vec2(1,1)); }
    /// Magnetic field at a point, with interaction matrix
    virtual vec2 fieldAtWithTransform2(const vec3& v, const Matrix<2,3,double>& M) const { return fieldAtWithTransform2(v,M,vec2(0,0),vec2(1,1)); }
    /// Magnetic field at a point, with interaction matrix
    virtual vec3 fieldAtWithTransform3(const vec3& v, const Matrix<3,3,double>& M) const { return fieldAtWithTransform3(v,M,vec2(0,0),vec2(1,1)); }
                
    /// Display field contributions over grid to given point
    void displayContribGrid(const vec3& v, unsigned int nx = 7, unsigned int ny = 7) const;
    /// draw a curve on the surface
    virtual void visualize_line(vec2 s, vec2 e) const;
    /// draw lines around a rectangular region
    virtual void visualize_region(vec2 ll, vec2 ur) const;
    
    SurfaceGeometry* mySurface; ///< surface on which source is defined

    unsigned int vis_n1;        ///< visualization gridding, z
    unsigned int vis_n2;        ///< visualization gridding, phi
};


#endif
