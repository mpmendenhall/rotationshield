#include "DipolePlane.hh"

bool DipolePlane::doNearfield = true;
bool DipolePlane::alwaysDipole = false;
unsigned int DipolePlane::drawSurfaces = 6;
unsigned int DipolePlane::nSurfaces = 6;

void DipolePlane::genSurfaces() {
	
	surfaces[0] = new PlaneSource(Plane(p.o+p.sn*d*0.5,p.dx,p.dz),murel);
	surfaces[1] = new PlaneSource(Plane(p.o-p.sn*d*0.5,-p.dx,-p.dz),murel);
	surfaces[2] = new PlaneSource(Plane(p.o+p.dx*0.5,-p.dz,-p.sn*d),murel);	// takes dx current << z state
	surfaces[3] = new PlaneSource(Plane(p.o-p.dx*0.5,p.dz,p.sn*d),murel);	// takes dx current << z state
	surfaces[4] = new PlaneSource(Plane(p.o+p.dz*0.5,-p.dx,-p.sn*d),murel);	// takes dz current << x state
	surfaces[5] = new PlaneSource(Plane(p.o-p.dz*0.5,p.dx,p.sn*d),murel);	// takes dz current << x state
	for(unsigned int i=0; i<6; i++)
		surfaces[i]->retain();
	
	//kfactor calculation
	
	/*
	 // kfactor estimation to first order in d
	 mdouble deltaz = (2*p.wx*p.wx+0*p.wz*p.wz)/(PI*p.wx*p.wz*sqrt(p.wx*p.wx+p.wz*p.wz))*d;
	 mdouble deltax = (0*p.wx*p.wx+2*p.wz*p.wz)/(PI*p.wx*p.wz*sqrt(p.wx*p.wx+p.wz*p.wz))*d;
	 kfactor = vec2( (mu-1)/(1+(mu-1)*deltax), (mu-1)/(1+(mu-1)*deltaz) )*2*d;
	 std::cout << kfactor;
	 */
	
	// direct kfactor determination from actual fields
	setState(vec3(1,0,0));
	vec2 bOutX = p.localProjection( fieldAt(p.o+p.sn*(0.500001*d)) );
	vec2 bInX = p.localProjection( fieldAt(p.o+p.sn*(0.499999*d)) );
	setState(vec3(0,1,0));
	vec2 bOutZ = p.localProjection( fieldAt(p.o+p.sn*(0.500001*d)) );
	vec2 bInZ = p.localProjection( fieldAt(p.o+p.sn*(0.499999*d)) );
	kfactor = vec2( (murel-1)/(bInX[0]-murel*bOutX[0]), (murel-1)/(bInZ[1]-murel*bOutZ[1]) );
	
	setState(state);	
}

Vec<6,mdouble> integral_components(mdouble x0, mdouble y0, mdouble z) {
	
	mdouble xx = x0*x0;	// x^2
	mdouble yy = y0*y0;	// y^2
	mdouble xy = x0*y0;	// x*y
	mdouble zz = z*z;	// z^2
	mdouble ir2xz = 1.0/(xx+zz);
	mdouble ir2yz = 1.0/(yy+zz);
	mdouble irxyz = 1.0/sqrt(xx + yy + zz);
	
	Vec<6,mdouble> v;
	
	v[0] = -xy*irxyz*ir2xz;				// xhat from mx
	v[1] = v[3] = irxyz;				// yhat from mx, xhat from my
	v[2] = z * ir2xz * (1.0-y0*irxyz); 	// zhat from mx
	v[4] = -xy*irxyz*ir2yz;				// yhat from my
	v[5] = z * ir2yz * (1.0-x0*irxyz);	// zhat from my
	
	return v/(4*PI);
}

mat3 DipolePlane::fieldAtComponentsNearfield(vec3 p0) const {
	/*
	 xsourced = zsourced = vec3();
	 
	 vec3 xs,zs;
	 
	 surfaces[0]->fieldAtComponents(p0,xs,zs);
	 xsourced -= zs;
	 zsourced += xs;
	 surfaces[1]->fieldAtComponents(p0,xs,zs);
	 xsourced -= zs;
	 zsourced += xs;
	 
	 if(nSurfaces == 6) {
	 surfaces[2]->fieldAtComponents(p0,xs,zs);
	 zsourced += zs;
	 surfaces[3]->fieldAtComponents(p0,xs,zs);
	 zsourced += zs;
	 
	 surfaces[4]->fieldAtComponents(p0,xs,zs);
	 xsourced -= zs;
	 surfaces[5]->fieldAtComponents(p0,xs,zs);
	 xsourced -= zs;
	 }
	 
	 xsourced /= 2*d;
	 zsourced /= 2*d;
	 */
	
	assert(false); //TODO
}


mat3 DipolePlane::fieldAtComponents(vec3 p0) const {
	
	assert(false); //TODO
	
	/*
	 
	 vec3 v = p.o - p0;				// vector from p0 to plane origin
	 mdouble z = v.dot(p.sn);		// distance from pt to plane
	 mdouble x0 = v.dot(p.dx)/p.wx;	// plane x center in rotated coordinates 
	 mdouble y0 = v.dot(p.dz)/p.wz;	// plane y center in rotated coordinates 
	 
	 if( ( doNearfield || ( fabs(z) < 50*d && fabs(x0) < (0.5*p.wx+50*d) && fabs(y0) < (0.5*p.wz+50*d) ) ) && !alwaysDipole ) {
	 
	 fieldAtComponentsNearfield(p0,xsourced,zsourced);
	 
	 } else {
	 
	 // farfield dipole plane approximation
	 
	 mdouble x1 = x0-0.5*p.wx;		// plane lower x bound in rotated coordinates
	 mdouble x2 = x0+0.5*p.wx;		// plane upper x bound in rotated coordinates
	 mdouble y1 = y0-0.5*p.wz;		// plane lower y bound in rotated coordinates
	 mdouble y2 = y0+0.5*p.wz;		// plane upper y bound in rotated coordinates
	 
	 Vec<6,mdouble> I = ( integral_components(x2,y2,z) - integral_components(x1,y2,z)
	 - integral_components(x2,y1,z) + integral_components(x1,y1,z) );
	 
	 xsourced = p.localToGlobal(vec3(I[0],I[1],I[2]));
	 zsourced = p.localToGlobal(vec3(I[3],I[4],I[5]));
	 
	 }
	 */
	
}
