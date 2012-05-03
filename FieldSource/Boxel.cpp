#include "Boxel.hh"
#include "Planesource.hh"

unsigned int Boxel::drawSurfaces = 2;

Boxel::Boxel(annulusSpec a, mdouble mu, mdouble thickness): CompoundPlane(Plane(a)), murel(mu), d(thickness) {
	
	annulusSpec a1 = a.translated(vec2(0,0.5*d)).flipped();
	Plane p1 = Plane(a1);
	annulusSpec a2 =  a.translated(vec2(0,-0.5*d));
	Plane p2 = Plane(a2);
	append(new PlaneSource(p1,murel));
	append(new PlaneSource(p2,murel));
	
	Plane p3 = Plane( annulusSpec(a1.end,a2.start,a.dTheta,a.theta0) );
	Plane p4 = Plane( annulusSpec(a2.end,a1.start,a.dTheta,a.theta0) );
	append(new PlaneSource(p3,murel));
	append(new PlaneSource(p4,murel));	
	
	Plane p5 = Plane(p.o+p.dx*0.5,(p1.o+p1.dx*0.5-p2.o-p2.dx*0.5),p.dz);
	append(new PlaneSource(p5,murel));
	
	Plane p6 = Plane(p.o-p.dx*0.5,(p1.o-p1.dx*0.5-p2.o+p2.dx*0.5),-p.dz);
	append(new PlaneSource(p6,murel));
	

}