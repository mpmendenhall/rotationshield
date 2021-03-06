/* 
 * FieldAdaptiveSurface.cpp, part of the RotationShield program
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

#include "FieldAdaptiveSurface.hh"
#include <iostream>

FieldAdaptiveSurface::FieldAdaptiveSurface(const DVFunc1<2,double>& f): l_remap(100), F(f) {
    setConstantSpacing();
    period = F.period;
}

vec2 FieldAdaptiveSurface::operator()(double x) const {
    if(F.period) x=wrap(x);
    return F(l_remap(x));
}

vec2 FieldAdaptiveSurface::deriv(double x, bool normalized) const {
    if(F.period) x=wrap(x);
    if(normalized) return F.deriv(l_remap(x),normalized);
    return F.deriv(l_remap(x)) * l_remap.deriv(x);
}

void FieldAdaptiveSurface::setConstantSpacing() {
    for(unsigned int i=0; i<l_remap.NX; i++)
        l_remap.set(i,double(i)/double(l_remap.NX-1));
}

void FieldAdaptiveSurface::symmetry_test() const {
    printf("Underlying function period: %g\n",F.period);
    for(int i=0; i<21; i++) {
        double l = float(i)/20.;
        double v = l_remap(l);
        printf("l=%g:\tv=%g\t%g\n",l,v,1-v);
    }
}

void FieldAdaptiveSurface::optimizeSpacing(const FieldEstimator2D& fes, double pfixed, bool useDeriv) {
    
    const unsigned int npts = 500;
    
    // cumulative parameter across length
    double fstr[npts];
    fstr[0]=0;
    for(unsigned int i=1; i < npts; i++) {
        double l = float(i-0.5)/float(npts-1.);
        vec2 zr = F(l);
        if(useDeriv) {
            // bunch up at high field derivatives along profile curve
            vec2 dv = F.deriv(l,true);
            fstr[i] = fstr[i-1] + fabs(fes.derivAt(zr,dv*0.0002).dot(dv));
        } else {
            // bunch up at high fields
            fstr[i] = fstr[i-1] + fes.estimateAt(zr).mag();
        }
    }
    
    // add "slope" to cumulative curve for fixed partitions
    pfixed = pfixed<=1e-3 ? 1e-3 : (pfixed>=0.999 ? 0.999:pfixed);
    double slope = (fstr[npts-1]+1e-8) * pfixed/(1.-pfixed) / float(npts);
    //printf("Optimized surface z spacing: to accumulated x=%g, adding slope %g\n",fstr[npts-1],slope*npts);
    for(unsigned int i=0; i < npts; i++) {
        fstr[i] += i*slope;
    }
    
    // convert to cumulative curve;
    for(unsigned int i=0; i<npts; i++)
        fstr[i] *= (l_remap.NX-1)/fstr[npts-1];
    
    // fill interpolator from inverse
    // with interpolation to determine dividing lines, relative coordinates 0 to 1
    l_remap.set(0,0.);
    unsigned int n=1;
    int i=1;
    while(n < l_remap.NX-1) {
        if(fstr[i] >= n) {
            double l = (double(i)-(fstr[i]-n)/(fstr[i]-fstr[i-1]))/(npts-1.); // interpolated partition point
            l_remap.set(n++, l);
            continue;
        }
        i++;
    }
    l_remap.set(l_remap.NX-1,1.);
}
