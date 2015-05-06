/* 
 * MixedSource.cpp, part of the RotationShield program
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

#include "MixedSource.hh"
#include "ScanRange.hh"

void MixedSource::loadSourcesFile(FILE* f, double scale = 1.0) {
    float x1,y1,z1,x2,y2,z2,j;
    while(fscanf(f,"LINE %f %f %f %f %f %f %f ",&x1,&y1,&z1,&x2,&y2,&z2,&j) == 7)
        addsource(new LineSource(vec3(x1,y1,z1),vec3(x2,y2,z2),j*scale));
}

vec3 MixedSource::fieldAt(const vec3& v) const {
    vec3 b(0,0,0);
    for(auto it = sources.begin(); it != sources.end(); it++) b += (*it)->fieldAt(v);
    return b;    
}

vec3 MixedSource::fieldOverLine(Line l) const {
    vec3 b = vec3();
    for(auto it = sources.begin(); it != sources.end(); it++) b += (*it)->fieldOverLine(l);
    return b;    
}

vec3 MixedSource::fieldOverPlane(Plane p) const {
    vec3 b = vec3();
    for(auto it = sources.begin(); it != sources.end(); it++) b += (*it)->fieldOverPlane(p);
    return b;    
}

void MixedSource::clear() {
    for(auto it = sources.begin(); it != sources.end(); it++) (*it)->release();
    sources.clear();
}

void MixedSource::arc(vec3 start, vec3 end, double j, int nsegs) {
    double r = sqrt(start[0]*start[0]+start[1]*start[1]);
    double z = start[2];
    double th0 = atan2(start[1],start[0]);
    double th1 = atan2(end[1],end[0]);
    if(th1-th0>M_PI) th0 += 2*M_PI;
    if(th0-th1>M_PI) th1 += 2*M_PI;
    if(th0 > th1)
    {
        double tmp = th1; th1 = th0; th0 = tmp;
        j = -j;
    }
    int nsteps = int((th1-th0)/(2.0*M_PI/double(nsegs))+1);
    double th = th0;
    vec3 v1,v2; v1[2] = v2[2] = z;
    v1[0] = r*cos(th); v1[1] = r*sin(th);
    for(int i=1; i<nsteps+1; i++)
    {
        th = th0+double(i)/double(nsteps)*(th1-th0);
        v2[0] = r*cos(th); v2[1] = r*sin(th);
        addsource(new LineSource(v1,v2,j));
        v1 = v2;
    }
}

void MixedSource::loop(double z, double r, double j, int nsides) {
    vec3 v0,v1;
    v0[2] = v1[2] = z;
    v0[0] = r; v0[1] = 0;
    ScanRange sr(0,2*M_PI,nsides);
    double th = sr.next();
    v0[0] = r*cos(th); v0[1] = r*sin(th);
    for(th = sr.next(); sr.goOn(); th = sr.next()) {
        v1[0] = r*cos(th); v1[1] = r*sin(th);
        addsource(new LineSource(v0,v1,j));
        v0 = v1;
    }
}

void MixedSource::addsource(const FieldSource* fs) {
    fs->retain();
    if(!sources.size()) mySymmetry = fs->getSymmetry();
    else mySymmetry += fs->getSymmetry();
    sources.push_back(fs);
}

void MixedSource::addsources(const MixedSource* MS) {
    for(auto it = MS->sources.begin(); it != MS->sources.end(); it++) addsource(*it);
}
