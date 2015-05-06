/* 
 * HoleDipolePerturbation.cpp, part of the RotationShield program
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

#include "HoleDipolePerturbation.hh"
#include "SurfaceCurrentRS.hh"

mvec HoleDipolePerturbation::getReactionTo(ReactiveSet* R, unsigned int phi) {
    assert(R && phi==0);
        
    BField_Protocol::BFP->x = x + mySurface.snorm(surfacePos,true)*dh;
    Matrix<2,3,double> RM2 =  Matrix<2,3,double>::identity() * mySurface.rotToLocal(surfacePos);
    BField_Protocol::BFP->M2 = &RM2;
    BField_Protocol::BFP->M2B = vec2(0,0);
    BField_Protocol::BFP->M3 = NULL;
    BField_Protocol::BFP->caller = RS_UID;
    if(!R->queryInteraction(BField_Protocol::BFP)) { assert(false); return mvec(); }
    
    return mvec(BField_Protocol::BFP->M2B);
}

bool HoleDipolePerturbation::queryInteraction(void* ip) {
    
    if(ip != BField_Protocol::BFP) return false;
    
    if(BField_Protocol::BFP->caller == RS_UID || hide_ixn.count(BField_Protocol::BFP->caller)) {
        return true;
    }
    
    if(BField_Protocol::BFP->M2) {
        BField_Protocol::BFP->M2B += fieldAtWithTransform2(BField_Protocol::BFP->x, *BField_Protocol::BFP->M2);
    } else if(BField_Protocol::BFP->M3) {
        BField_Protocol::BFP->B += fieldAtWithTransform3(BField_Protocol::BFP->x, *BField_Protocol::BFP->M3);
    } else {
        BField_Protocol::BFP->B += fieldAt(BField_Protocol::BFP->x);
    }
    return true;
}

void HoleDipolePerturbation::_setDF(unsigned int, double) {
    vec3 d0 = mySurface.deriv(surfacePos,0,true);
    vec3 d1 = mySurface.deriv(surfacePos,1,true);
    
    m = (d0*finalState[0] + d1*finalState[1]) * 8. * a*a*a / 3.;
}

void HoleDipolePerturbation::_visualize() const {
    DipoleSource::_visualize();
    vsr::setColor(0,0,1);
    vsr::dot(x + mySurface.snorm(surfacePos,true)*dh);
    
    vec3 d0 = mySurface.deriv(surfacePos,0,true);
    vec3 d1 = mySurface.deriv(surfacePos,1,true);
    vsr::setColor(1,0,0);
    vsr::startLines();
    for(int i=0; i<=21; i++) {
        double th = 2*M_PI*i/21.;
        vsr::vertex( x + d0*(a*cos(th)) + d1*(a*sin(th)) );
    }
    vsr::endLines();
}
