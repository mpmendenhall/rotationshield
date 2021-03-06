/* 
 * Geometry.cpp, part of the RotationShield program
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

#include "Geometry.hh"
#include <algorithm>

void Line::_visualize() const { vsr::line(s,e); }

void Line::visualizeDirected(float j) const { 
    vsr::startLines();
    for(float i=0; i<8; i++) {
        float p = i/7.0;
        float q = 1.0-p;
        if(p<0.5)
            vsr::setColor(2*p,0,1.0,1.0);
        else
            vsr::setColor(1.0,0,2*q,1.0);
        if(j>0)
            vsr::vertex(s*q+e*p);
        else
            vsr::vertex(s*p+e*q);
    }
    vsr::endLines();
}

Plane::Plane(annulusSpec a) {
    vec2 c = (a.start+a.end)*0.5;
    vec2 d = a.end-a.start;
    o = zrotate(vec3(0,c[1],c[0]),a.theta0);
    dx = zrotate(vec3(2*c[1]*tan(a.dTheta*0.5),0,0),a.theta0);
    dz = zrotate(vec3(0,d[1],d[0]),a.theta0);
    wx = dx.mag();
    wz = dz.mag();
    area = wx*wz;
    sn = cross(dx,dz)/area;
}

const mat3 Plane::projectionMatrix() const {
    mat3 M;
    for(unsigned int i=0; i<3; i++) {
        M(0,i) = dx[i]/wx;
        M(1,i) = dz[i]/wz;
        M(2,i) = sn[i];
    }
    return M;
}


void Plane::_visualize() const { vsr::plane(o,dx,dz); }

void Plane::visualizeCoordinates(float scale) const {
    float l = std::min(wx,wz)*0.5*scale;
    vsr::line(o, o+dx*l/wx);
    vsr::line(o, o+dz*l/wz);
    vsr::line(o, o+sn*l);
}

void Plane::visualizeVector(vec3 v) const {
    v = v.normalized()*(0.5*std::min(wx,wz));
    vsr::line(o,o+dx*v[0]/wx+dz*v[1]/wz+sn*v[2]);
}

Plane Plane::readFromFile(std::istream& s) { 
    vec3 origin = vec3::readBinary(s);
    vec3 dxx = vec3::readBinary(s);
    vec3 dzz = vec3::readBinary(s);
    return Plane(origin,dxx,dzz);
}


vec3 Plane::zrotate(const vec3& v, double th) {
    double cs = cos(th);
    double sn = sin(th);
    return vec3( cs*v[0]-sn*v[1], sn*v[0]+cs*v[1], v[2] );
}

vector<Plane> Plane::subdivide(unsigned int nX, unsigned int nZ) const {
    assert(nX>0 && nZ>0);
    vector<Plane> v = vector<Plane>();
    for(unsigned int x=0; x<nX; x++)
        for(unsigned int z=0; z<nZ; z++)
            v.push_back( Plane(o - dx*(0.5 - (double(x)+0.5)/nX) - dz*(0.5 - (double(z)+0.5)/nZ), dx/nX, dz/nZ) );
    return v;
}
