#include "PlaneSource.hh"
#include <math.h>

mmat PlaneSource::fieldAtComponents(vec3 p0) const {
	
	vec3 v = p0 - p.o;
	mdouble a = v.dot(p.sn);
	mdouble aa = a*a;
	mdouble x0 = v.dot(p.dx)/p.wx;
	mdouble x1 = x0-0.5*p.wx;
	mdouble x2 = x0+0.5*p.wx;
	mdouble y0 = v.dot(p.dz)/p.wz;
	mdouble y1 = y0-0.5*p.wz;
	mdouble y2 = y0+0.5*p.wz;
	
	// perpendicular component
	mdouble c11 = sqrt(aa+x1*x1+y1*y1);
	mdouble c12 = sqrt(aa+x1*x1+y2*y2);
	mdouble c21 = sqrt(aa+x2*x2+y1*y1);
	mdouble c22 = sqrt(aa+x2*x2+y2*y2);
	mdouble perp1 = -1.0/(8*M_PI)*log( (x1+c12)*(x2-c22)*(x1-c11)*(x2+c21) /
									((x1-c12)*(x2+c22)*(x1+c11)*(x2-c21)) );
	mdouble perp2 = 1.0/(8*M_PI)*log( (y2+c12)*(y2-c22)*(y1-c11)*(y1+c21) /
								   ((y2-c12)*(y2+c22)*(y1+c11)*(y1-c21)) );
	
	// parallel component
	mdouble z11 = x1*y1/(a*c11);
	mdouble z12 = x1*y2/(a*c12);
	mdouble z21 = x2*y1/(a*c21);
	mdouble z22 = x2*y2/(a*c22);
	mdouble par = 1.0/(4*M_PI)*(atan(z22)-atan(z21)+atan(z11)-atan(z12));
	
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
