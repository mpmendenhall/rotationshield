/* 
 * PlanarElement.cpp, part of the RotationShield program
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

#include "PlanarElement.hh"
#include <algorithm>

void PlanarElement::_visualize() const {

    double logmax = 3.0;
    float smag = state.mag2();
    if(smag)
        smag = std::max(0.,std::min(1.0,0.1*(log(smag)+10-logmax)));
    
    vec3 svec = vec3();
    for(unsigned int i=0; i<std::min((unsigned int)3, state.size()); i++)
        svec[i] = state[i];
    
    vec3 hsv = vec3( atan2(svec[0],svec[1]), 1.0, 1.0 );
    vec3 rgb = hsv2rgb(hsv);
    vsr::setColor(rgb[0], rgb[1] , rgb[2], 1.0*smag);
    p._visualize();
    vsr::setColor(1, .7, .7, 1.0);
    p.visualizeCoordinates(0.2);
    vsr::setColor(.7, .7, 1, 1.0);
    p.visualizeCoordinates(-0.2);
    vsr::setColor(1,1,1,1);
    
    p.visualizeVector(svec);
}
