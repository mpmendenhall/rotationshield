/* 
 * InterpolatingRS.cpp, part of the RotationShield program
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

#include "InterpolatingRS.hh"
#include <cassert>
#include <stdio.h>

void InterpolatingRS2D::clear_data() {
    for(unsigned int i=0; i<G.size(); i++)
        delete G[i];
    G.clear();
}

void InterpolatingRS2D::make_grids(unsigned int nz, unsigned int ndf) {
    clear_data();
    nDFi = ndf;
    nZ = nz;
    for(unsigned int i=0; i<nDFi; i++) {
        G.push_back(new BicubicGrid(nZ,nPhi));
        if(nPhi>1) {
            G.back()->bc[1] = IB_CYCLIC;
            G.back()->setUserRange(0,1,true,0.5);
            G.back()->setUserRange(0,1,false,0.5);
        }
    }
}

void InterpolatingRS2D::DF_address(unsigned int DF, unsigned int& p, unsigned int& z, unsigned int& d) const {
    p = DF%nPhi;
    z = (DF/nPhi)%nZ;
    d = DF/(nPhi*nZ);
}

void InterpolatingRS2D::_setDF(unsigned int DF, double v) {
    assert(DF < nDF());
    assert(v==v && fabs(v)<1e12);
    unsigned int p,z,d;
    DF_address(DF,p,z,d);
    assert(d<G.size() && p<nPhi && z<nZ);
    G[d]->set(z,p,v);
}

mvec InterpolatingRS2D::interpl_DF(vec2 l) const {
    mvec v(nDFi);
    for(unsigned int i=0; i<nDFi; i++)
        v[i] = (*G[i])(l[0],l[1]);
    return v;
}

void InterpolatingRS2D::printData() const {
    for(unsigned int i=0; i<nDFi; i++) {
        printf("Data table for layer %i:\n",i);
        G[i]->printData();
    }
}
