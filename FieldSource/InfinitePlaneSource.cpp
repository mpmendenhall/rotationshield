#include "InfinitePlaneSource.hh"

mmat InfinitePlaneSource::fieldAtComponents(vec3 p0) const {
	
	vec3 v = p0 - p.o;
	mdouble a = v.dot(p.sn);
	mdouble aa = a*a;
	mdouble x0 = v.dot(p.dx)/p.wx;
	mdouble x1 = x0-0.5*p.wx;
	mdouble x2 = x0+0.5*p.wx;
	
	// parallel components
	mdouble par = (atan(x2/a)-atan(x1/a))/(2*PI);
	
	// perpendicular component
	mdouble perp = (log((x2*x2+aa)/(x1*x1+aa)))/(4*PI);
	
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
