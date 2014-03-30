/* 
 * PlaneSource.cpp, part of the RotationShield program
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

#include "PlaneSource.hh"
#include <math.h>

mmat PlaneSource::fieldAtComponents(vec3 p0) const {
	
	vec3 v = p0 - p.o;
	double a = v.dot(p.sn);
	double aa = a*a;
	double x0 = v.dot(p.dx)/p.wx;
	double x1 = x0-0.5*p.wx;
	double x2 = x0+0.5*p.wx;
	double y0 = v.dot(p.dz)/p.wz;
	double y1 = y0-0.5*p.wz;
	double y2 = y0+0.5*p.wz;
	
	// perpendicular component
	double c11 = sqrt(aa+x1*x1+y1*y1);
	double c12 = sqrt(aa+x1*x1+y2*y2);
	double c21 = sqrt(aa+x2*x2+y1*y1);
	double c22 = sqrt(aa+x2*x2+y2*y2);
	double perp1 = -1.0/(8*M_PI)*log( (x1+c12)*(x2-c22)*(x1-c11)*(x2+c21) /
									((x1-c12)*(x2+c22)*(x1+c11)*(x2-c21)) );
	double perp2 = 1.0/(8*M_PI)*log( (y2+c12)*(y2-c22)*(y1-c11)*(y1+c21) /
								   ((y2-c12)*(y2+c22)*(y1+c11)*(y1-c21)) );
	
	// parallel component
	double z11 = x1*y1/(a*c11);
	double z12 = x1*y2/(a*c12);
	double z21 = x2*y1/(a*c21);
	double z22 = x2*y2/(a*c22);
	double par = 1.0/(4*M_PI)*(atan(z22)-atan(z21)+atan(z11)-atan(z12));
	
	// filter out NANs
	if(!(perp1 == perp1)) perp1 = 0;
	if(!(perp2 == perp2)) perp2 = 0;
	if(!(par == par)) par = 0;
	
	vec3 xsourced = p.sn*perp1 + p.dz*par/p.wz;
	vec3 zsourced = -(p.sn*perp2 - p.dx*par/p.wx);
	
	mmat f = mmat(3,2);
	for(unsigned int i=0; i<3; i++) {
		f(i,0) = xsourced[i];
		f(i,1) = zsourced[i];
	}
	
	return f;
	
}
