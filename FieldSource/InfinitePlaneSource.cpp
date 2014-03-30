/* 
 * InfinitePlaneSource.cpp, part of the RotationShield program
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

#include "InfinitePlaneSource.hh"
#include <math.h>

mmat InfinitePlaneSource::fieldAtComponents(vec3 p0) const {
	
	vec3 v = p0 - p.o;
	double a = v.dot(p.sn);
	double aa = a*a;
	double x0 = v.dot(p.dx)/p.wx;
	double x1 = x0-0.5*p.wx;
	double x2 = x0+0.5*p.wx;
	
	// parallel components
	double par = (atan(x2/a)-atan(x1/a))/(2*M_PI);
	
	// perpendicular component
	double perp = (log((x2*x2+aa)/(x1*x1+aa)))/(4*M_PI);
	
	// filter out NANs
	if(!(perp == perp)) perp = 0;
	
	vec3 xsourced = vec3(); //p.dz*par/p.wz;
	vec3 zsourced = p.sn*perp - p.dx*par/p.wx;
	
	
	mmat f = mmat(3,2);
	for(unsigned int i=0; i<3; i++) {
		f(i,0) = -xsourced[i];
		f(i,1) = -1.0218*zsourced[i];
	}
	
	return f;
	
}
