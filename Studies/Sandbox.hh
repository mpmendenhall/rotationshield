#ifndef SANDBOX_HH
#define SANDBOX_HH 1

#include "SymmetricSolver.hh"
#include "GenericSolver.hh"
#include "FieldSource.hh"
#include "analysis.hh"
#include "SurfacelCyl.hh"
#include "CosThetaBuilder.hh"
#include "ScanRange.hh"
#include <iostream>
#include "Planesource.hh"
#include "InfinitePlaneSource.hh"
#include "UniformField.hh"
//#include "CompoundElement.hh"
//#include "DipolePlane.hh"
#include "Boxel.hh"
//#include "GridBoxel.hh"
#include "tests.hh"

//octant of sample volume
vec3 sampleOctantLL(0.05,-0.051,0);
vec3 sampleOctantUR(0.126,0.0,0.25);
//entire sample volume
vec3 sampleLL(0.05,-0.051,-0.25);
vec3 sampleUR(0.126,0.051,0.25);


void comparison_shield(std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	
	CosThetaBuilder(17, 0.61, 3.92).buildCoil(*MxS);
	
	
	FieldEstimator2D fe;
	fe.addsource(vec2(-3.92/2,0.61),1.0);
	fe.addsource(vec2(3.92/2,0.61),1.0);
	
	
	SurfacelCyl SB = SurfacelCyl(128);
	//SB.makeOptCyl(20, 30, .6223, -3.9624/2, 3.9624/2, new DipolePlane(Plane(),100000.0,0.001), &fe);
	
	SymmetricSolver* s = new SymmetricSolver(&SB);
	s->retain();
	s->visualize();
	s->calculateIncident(MxS);
	s->visualize();
	s->solve();
	s->calculateResult();
	
	MxS->addsource(s);
	
	MxS->visualize();
	
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	FA.survey(vec3(-.15,-.06,-.25),vec3(.15,.06,.25),9,9,9,fitf,gridf);
	
	vsr::pause();
	
	s->release();
	MxS->release();
}


void shieldingFactorTest(std::ostream& outf, mdouble mu) {
	
	SurfacelCyl SB = SurfacelCyl(128);
	//SB.makeOptCyl(30, 0, .6223, -3.9624/2, 3.9624/2, new DipolePlane(Plane(),mu,0.001) );

	SymmetricSolver* s = new SymmetricSolver(&SB);
	s->retain();
	s->visualize();
	
	s->solve();
	

	
	s->calculateIncident(new UniformField(vec3(1,0,0)));
	s->calculateResult();
	s->visualize();
	vsr::pause();
	
	MixedSource* MxS = new MixedSource();
	MxS->addsource(s);
	MxS->addsource(new UniformField(vec3(1,0,0)));
	
	outf << mu << "\t" << MxS->fieldAt(vec3(0,0,0))[0] << "\t";
	
	s->calculateIncident(new UniformField(vec3(0,0,1)));
	s->calculateResult();
	s->visualize();
	
	MixedSource* MxSb = new MixedSource();
	MxSb->addsource(s);
	MxSb->addsource(new UniformField(vec3(0,0,1)));
	outf << MxSb->fieldAt(vec3(0,0,0))[2] << "\n";
	
	vsr::pause();
	
	s->release();
}


/// Test DipolePlane properties
void dipoleTest() {
	
	//Doublewall* w = new Doublewall(annulusSpec(vec2(0,0.5),vec2(0.5,0.5),1.0,0),1e5,0.01);
	//GridBoxel* w = new GridBoxel(annulusSpec(vec2(0,0.5),vec2(0.5,0.5),1.0,0),1e4,0.1,8,5,4);
	Boxel* w = new Boxel(annulusSpec(vec2(0,0.5),vec2(0.5,0.5),1.0,0),1e8,0.0001);
	//PlanarElement* w = 0;
	
	//Boxel* w = new Boxel(annulusSpec(vec2(0,0.5),vec2(0.5,0.5),1.0,0),1e9,0.001);
	w->retain();
	
	unsigned int nEls = 128;
	unsigned int nZ = 20;
	SurfacelCyl SB = SurfacelCyl(nEls);
	
	if(true) {
		SB.makeOptCyl(nZ, 0, .50, -1, 1, w);
	} else {
		FieldEstimator2D fe;
		fe.addsource(vec2(-2,0.40),1.0);
		fe.addsource(vec2(2,0.40),1.0);
		SB.makeOptCyl(11, 10, .55, -2, 2, w, &fe);
		SB.makeOptCyl(11, 10, .70, -2, 2, w, &fe);
	}
	
	InteractionSolver* G;
	
	if(true) {
		G = new SymmetricSolver(&SB);
	} else {
		G = new GenericSolver();
		for(unsigned int i=0; i<1; i++) {
			G->addSurfacel(SB.genElement(nZ/2, i));
		}
		
	}
	G->retain();
	G->solve();
	
	UniformField* f0 = new UniformField(vec3(1.0,0.0,0.0));
	f0->retain();
	
	G->calculateIncident(f0);
	G->calculateResult();
	
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	MxS->addsource(G);
	MxS->addsource(f0);
	
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	FA.visualizeSurvey(vec3(-1,-1,0), vec3(1,1,0), 31,31,1);
	//FA.visualizeSurvey(vec3(-0.5,0,0), vec3(0.5,0.75,0), 31,31,1);
	//FA.visualizeSurvey(vec3(0,0,-0.5), vec3(0,0.75,0.5), 1,31,31);


	vsr::pause();
}

void infiniTest() {
	
	InfinitePlaneSource* w = new InfinitePlaneSource(vec2(),vec2(),1e20);
	w->retain();
	
	unsigned int nEls = 64;
	unsigned int nZ = 1;
	SurfacelCyl SB = SurfacelCyl(nEls);
	
	if(false) {
		SB.makeOptCyl(nZ, 0, .50, 1, -1, w);
	} else {
		SB.makeOptCyl(nZ, 0, .45, -1, 1, w);
		SB.makeOptCyl(nZ, 0, 0.55, 1, -1, w);
	}
	
	InteractionSolver* G;
	
	if(true) {
		G = new SymmetricSolver(&SB);
	} else {
		G = new GenericSolver();
		for(unsigned int i=0; i<1; i++) {
			G->addSurfacel(SB.genElement(nZ/2, i));
		}
		
	}
	G->retain();
	G->solve();
	
	UniformField* f0 = new UniformField(vec3(1.0,0.0,0.0));
	f0->retain();
	
	G->calculateIncident(f0);
	G->calculateResult();
	
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	MxS->addsource(G);
	MxS->addsource(f0);
	
	//std::cout << MxS->fieldAt(vec3(0.1,0.499999999,0)) << std::endl;
	std::cout << MxS->fieldAt(vec3()) << std::endl;
	
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	FA.visualizeSurvey(vec3(-1,-1,0), vec3(1,1,0), 31,31,1);
	//FA.visualizeSurvey(vec3(-0.5,.45,0), vec3(0.5,0.55,0), 31,31,1);
	//FA.visualizeSurvey(vec3(0,0,-0.5), vec3(0,0.75,0.5), 1,31,31);
	
	
	vsr::pause();
}




#endif
