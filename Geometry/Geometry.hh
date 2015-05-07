/* 
 * Geometry.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
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

#ifndef GEOMETRY_HH
/// Make sure this header is only loaded once
#define GEOMETRY_HH

#include <cmath>
#include <algorithm>
#include "Visr.hh"
#include "Typedefs.hh"
#include <vector>
#include <stdio.h>

/// A (directed) line, represented by two endpoint vectors
class Line: public Visualizable {
public:
    /// Contructor
    Line(vec3 start, vec3 end): s(start), e(end), len((end-start).mag()), dv((end-start).normalized()) {}
    /// Default constructor (uninitialized endpoints)
    Line(): len(0) {}
    /// Get a point along the line, parametrized by x
    vec3 position(double x) const { return s + (e-s)*x; }
    ///print info to stdout
    void display() const { printf("Line: "); s.display(); printf(" to "); e.display("\n"); }
    ///draw to visualizer
    void _visualize() const;
    void visualizeDirected(float j = 1.0) const;
    
    const vec3 s;       ///< Starting point
    const vec3 e;       ///< Ending point
    const double len;   ///< Length of the line, \f$ |e-s| \f$
    vec3 dv;            ///< Normalized direction \f$ \frac{e-s}{|e-s|} \f$
};


/// specification for a segment of an annulus around the z axis
class annulusSpec {
public:
    /// constructor
    annulusSpec(vec2 s = vec2(), vec2 e = vec2(), double dth = 0, double th0 = 0): start(s), end(e), theta0(th0), dTheta(dth) {}
    
    /// length along z
    double dl() const { return (end-start).mag(); }
    /// central length along phi
    double dr() const { return (end[1]+start[1])*0.5*dTheta; }
    
    /// annulusSpec with translated start and endpoints
    annulusSpec translated(vec2 v) const { annulusSpec a = *this; a.start += v; a.end += v; return a; }
    
    /// annulusSpec with start/end swapped
    annulusSpec flipped() const { annulusSpec a = *this; std::swap(a.start,a.end); return a; }
    
    /// print to stdout
    void display() const { std::cout << "annulusSpec "<< start << "-" << end << " " << dTheta << "@" << theta0 << std::endl; }
    
    /// split into grid of sub-sections
    vector<annulusSpec> subdivide(unsigned int nZ, unsigned int nTheta) const {
        assert(nZ>0 && nTheta>0);
        vector<annulusSpec> v = vector<annulusSpec>();
        for(unsigned int z=0; z<nZ; z++) {
            double lz0 = 1.0-double(z)/double(nZ);
            double lz1 = 1.0-double(z+1)/double(nZ);
            for(unsigned int t=0; t<nTheta; t++)
                v.push_back( annulusSpec(start*lz0+end*(1.0-lz0),start*lz1+end*(1.0-lz1), dTheta/double(nTheta),
                                         theta0 - (0.5 - (double(t)+0.5)/nTheta)*dTheta ) );
        }
        return v;
    }
    
    vec2 start;         ///< starting (z,r)
    vec2 end;           ///< ending (z,r)
    double theta0;      ///< center angle
    double dTheta;      ///< subtended angle
};

/// An (oriented) plane, represented by an origin and two width/direction vectors
class Plane: public Visualizable {
public:
    /// Constructor
    Plane(vec3 origin = vec3(), vec3 dxx = vec3(1,0,0), vec3 dzz = vec3(0,0,1)): wx(dxx.mag()), wz(dzz.mag()), o(origin), dx(dxx), dz(dzz)  { area = wx*wz; sn = cross(dx,dz)/area; }
    /// Constructor from annulusSpec
    Plane(annulusSpec a);
    
    /// Parametrizes positions on the plane in a local coordinate system
    vec3 position(double x, double z) const { return o + dx*(x*0.5) + dz*(0.5*z); }
    /// Project a vec3 onto the plane's local coordinates
    vec2 localProjection(vec3 v) const { return vec2(v.dot(dx)/wx, v.dot(dz)/wz); }
    /// get projection matrix to local coordinates
    const mat3 projectionMatrix() const;
    /// Print info about the plane to stdout
    void display() const { printf("Plane (a=%g)\n\to:\t",(double)area); o.display(); printf("\tdx:\t"); dx.display(); printf("\tdz:\t"); dz.display(); }
    
    /// get subtended angle around z axis
    double dTheta() const { return 2*atan2(0.5*wx,vec2(o[0],o[1]).mag()); }
    ///draw to visualizer
    void _visualize() const;
    ///draw local coordinate axes
    void visualizeCoordinates(float scale = 1.0) const;
    /// visualize a vector in local coordinates
    void visualizeVector(vec3 v) const;
    /// write to binary output file
    void writeToFile(std::ostream& of) const { o.writeBinary(of); dx.writeBinary(of); dz.writeBinary(of); }
    /// read from binary output file
    static Plane readFromFile(std::istream& s);
    /// transform local plane-coordinates vector to global coordinates
    vec3 localToGlobal(const vec3& v) const { return dx*(v[0]/wx) + dz*(v[1]/wz) + sn*v[2]; }
    /// transform local plane-coordinates vector to global coordinates
    vec3 surfacePosition(const vec2& v) const { return o + dx*v[0] + dz*v[1]; }
    /// plane rotated around the z axis
    Plane zrotated(double th) const { return Plane(zrotate(o,th),zrotate(dx,th),zrotate(dz,th)); }
    /// Rotate a vector around the z axis
    static vec3 zrotate(const vec3& v, double th);
    /// split into grid of sub-sections
    vector<Plane> subdivide(unsigned int nX, unsigned int nZ) const;
    
    double area;        ///< area of the plane
    double wx;          ///< width in "x" (\f$\hat\phi\f$) direction
    double wz;          ///< width along z direction
    vec3 o;             ///< center of the plane
    vec3 dx;            ///< vector along the "x" edge of the plane
    vec3 dz;            ///< vector along the "z" edge of the plane
    vec3 sn;            ///< Surface normal \f$ \frac{dx \times dz}{|dx||dz|} \f$
    
};

class AnnularBoxel: public annulusSpec {
public:
    AnnularBoxel(annulusSpec a, double thickness): annulusSpec(a), d(thickness) {
        top = translated(vec2(0,0.5*d)).flipped();
        bottom = translated(vec2(0,-0.5*d));
        front = annulusSpec(top.end,bottom.start,dTheta,theta0);
        back = annulusSpec(bottom.end,top.start,dTheta,theta0);
        Plane p = Plane(*this);
        Plane p1 = Plane(top);
        Plane p2 = Plane(bottom);
        left = Plane(p.o+p.dx*0.5,(p1.o+p1.dx*0.5-p2.o-p2.dx*0.5),p.dz);
        right = Plane(p.o-p.dx*0.5,(p1.o-p1.dx*0.5-p2.o+p2.dx*0.5),-p.dz);
    }
    
    double d;
    annulusSpec top;
    annulusSpec bottom;
    annulusSpec front;
    annulusSpec back;
    Plane left;
    Plane right;
};



#endif
