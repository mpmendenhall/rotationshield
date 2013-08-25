#include "LineSource.hh"
#include "gsl/gsl_sf_ellint.h"
#include <math.h>

vec3 LineSource::fieldAt(const vec3& v) const
{
	//difference vectors to each endpoint
	vec3 dv1 = l.s - v;
	vec3 dv2 = l.e - v;
	
	//field calculation
	vec3 b = cross(dv1,l.dv);
	mdouble m2 = b.mag2();
	if(!m2) return b;
	b *= j*(dv2.dot(l.dv)/dv2.mag() - dv1.dot(l.dv)/dv1.mag())/(4.0*M_PI*m2);
	return b;	
}


mdouble LineSource::sF1(mdouble a2, mdouble b2, mdouble x)
{
	if(x<0) return -sF1(a2,b2,-x);
	mdouble phi = atan(x/sqrt(b2));
	mdouble k = sqrt((a2-b2)/a2); // = sqrt(m) in A&S notation
	mdouble a = gsl_sf_ellint_E(phi, k, GSL_PREC_DOUBLE)*sqrt(a2)/b2; // other modes: GSL_PREC_DOUBLE, GSL_PREC_APPROX
	//printf("SF1(%g,%g,%g) = %g\n",a2,b2,x,a);
	return a;
}
mdouble LineSource::sF2(mdouble a2, mdouble b2, mdouble x)
{
	if(x<0) return -sF2(a2,b2,-x);
	mdouble phi = asin(sqrt(x*x*(b2-a2)/(b2*(a2+x*x))));
	mdouble k = sqrt(b2/(b2-a2)); // = sqrt(m) in A&S notation
	return gsl_sf_ellint_E(phi, k, GSL_PREC_SINGLE)*a2*sqrt(a2-b2); // other modes: GSL_PREC_DOUBLE, GSL_PREC_APPROX
}

/// Integrates \f$ \int_{y_1}^{y_2} \frac{\sqrt{a^2+b^2+y^2}}{a^2+y^2}dy \f$
mdouble sF3(mdouble a2, mdouble b2, mdouble y1, mdouble y2)
{
	mdouble c1 = sqrt((a2-b2)/b2); 
	return log((y2+sqrt(a2+y2*y2))/(y1+sqrt(a2+y1*y1))) +c1*( atan(c1*y2/sqrt(a2+y2*y2)) - atan(c1*y1/sqrt(a2+y1*y1)) );
}

vec3 LineSource::fieldOverLine(Line L) const
{
	mdouble relDirection = l.dv.dot(L.dv);
	if(fabs(1-relDirection) < 1e-8) {
		//     x3-----------x4 L
		// v -/ |- a
		//   0-------x2 l
		vec3 v = L.s - l.s;
		mdouble x2 = l.len;
		mdouble x3 = v.dot(l.dv);
		mdouble x4 = x3 + L.len*relDirection;
		mdouble aa = v.mag2() - x3*x3;
		if(aa == 0) { return vec3(); }
		mdouble b = relDirection*(sqrt(aa+x4*x4)+sqrt(aa+(x2-x3)*(x2-x3))-sqrt(aa+x3*x3)-sqrt(aa+(x2-x4)*(x2-x4)))/sqrt(aa);
		vec3 d = cross(l.dv,v);
		return d*(j*b/(4.0*M_PI*d.mag()*L.len));
	} else if(fabs(relDirection)<1e-8)
	{
		//         y2
		//     l   /                  . L
		//   x1---O-----x2    a-[ x1________x2 l
		//    \  /
		//   v \/ L
		//     y1
		
		vec3 v = L.s - l.s;
		mdouble x1 = -v.dot(l.dv);
		mdouble x2 = x1+l.len;
		mdouble y1 = -v.dot(L.dv);
		mdouble y2 = y1 - L.len;
		mdouble aa = v.mag2() - x1*x1 - y1*y1;
		vec3 z = cross(l.dv,L.dv);
		mdouble a = v.dot(z);
		
		// perpendicular component
		mdouble c11 = sqrt(aa+x1*x1+y1*y1);
		mdouble c12 = sqrt(aa+x1*x1+y2*y2);
		mdouble c21 = sqrt(aa+x2*x2+y1*y1);
		mdouble c22 = sqrt(aa+x2*x2+y2*y2);
		mdouble perp = j/(8*M_PI*L.len)*log( (x1+c12)*(x2-c22)*(x1-c11)*(x2+c21) /
										  ((x1-c12)*(x2+c22)*(x1+c11)*(x2-c21)) );
		
		// parallel component
		mdouble z11 = x1*y1/(a*c11);
		mdouble z12 = x1*y2/(a*c12);
		mdouble z21 = x2*y1/(a*c21);
		mdouble z22 = x2*y2/(a*c22);
		mdouble par = j/(4*M_PI*L.len)*(atan(z22)-atan(z21)+atan(z11)-atan(z12));
		
		return z*perp + L.dv*par;
	}
	return FieldSource::fieldOverLine(L);
}


/*
vec2 LineSource::surfaceField(Plane p) const {
	mdouble relDirection = l.dv.dot(p.dz)/p.wz;
	
	if(fabs(relDirection)>0.99999) {
		
		//       y1
		//       +============+
		//   x3 /===== p ====/ x4 
		//     +===========+ 
		//    y2  v -/ |- a
		//          0-------x2 l
		
		vec3 v = p.o - l.s;
		mdouble a = v.dot(p.sn);
		mdouble aa = a*a;
		mdouble x0 = v.dot(l.dv);
		mdouble x2 = l.len;
		mdouble x3 = x0 - 0.5*p.wz;
		mdouble x4 = x0 + 0.5*p.wz;
		mdouble y0 = v.dot(p.dx)/p.wx;
		mdouble y1 = y0 - 0.5*p.wx;
		mdouble y2 = y0 + 0.5*p.wx;
		vec2 r; r[1] = 0;
		r[0] = -relDirection*( sF3(aa+(x2-x3)*(x2-x3),aa,y1,y2) - sF3(aa+x3*x3,aa,y1,y2)
							  +sF3(aa+x4*x4,aa,y1,y2) - sF3(aa+(x2-x4)*(x2-x4),aa,y1,y2) )*fabs(a)*j/(4.0*M_PI*p.area);
		return r;
	}
}
 */
