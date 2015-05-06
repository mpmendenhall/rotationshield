/* 
 * FieldEstimator2D.cpp, part of the RotationShield program
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

#include "FieldEstimator2D.hh"

vec2 FieldEstimator2D::estimateAt(const vec2& v) const {
    vec2 est = vec2();
    vec2 r;
    double r0;
    for(unsigned int i=0; i<sources.size(); i++) {
        r = v-sources[i];
        r0 = r.mag2();
        est[0] += r[1]*(-currents[i]/r0);
        est[1] += r[0]*(currents[i]/r0);
    }
    return est;
}

vec2 FieldEstimator2D::derivAt(const vec2& v, const vec2& dx) const {
    return (estimateAt(v+dx*0.5)-estimateAt(v-dx*0.5))/dx.mag();
}

vec2 FieldEstimator2Dfrom3D::estimateAt(const vec2& v) const {
    vec3 F = myFS->fieldAt(vec3(v[1],0,v[0]));
    return vec2(F[2],F[0]);
}
