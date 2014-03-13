#include "tests.hh"
#include "SymmetricSolver.hh"
#include "GenericSolver.hh"
#include "FieldSource.hh"
#include "PlaneSource.hh"
#include "SurfacelCyl.hh"
#include "analysis.hh"
#include "CosThetaBuilder.hh"
#include "Integrator.hh"
#include "SurfaceCurrentSource.hh"
#include "SurfaceCurrentRS.hh"
#include "FieldAdaptiveSurface.hh"

bool compareResults(mdouble a, mdouble b, const char* label) {
	bool pass = true;
	if(label) printf("%s:\n",label);
	printf("\tTest Result:\t%.14e\n\t  Reference:\t%.14e\n\t\t   error:\t%.1f%%",double(a),double(b),double(100*(a-b)/a));
	if(fabs(100*(a-b)/a) > 5.0 ) { printf(" <<-------- ******* TEST FAILED ******"); pass = false; }
	printf("\n");
	fflush(stdout);
	return pass;
}

bool reference_sanity_check() {
	printf("Testing bare wire fields...\n");
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	
	CosThetaBuilder CB(17, 0.61, 3.92);
	CB.myCap[0] = CB.myCap[1] = CosThetaBuilder::CAP_LINE;
	CB.buildCoil(*MxS);
	MxS->visualize();
	
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan = vec3(0.15,0,0);
	vec3 yscan = vec3(0,0.06,0);
	vec3 zscan = vec3(0,0,0.25);
	
	bool pass = true;
	mdouble b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,8.53210605848347e-01,"coil origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,8.53217812021647e-01-8.53210605848347e-01,"coil +x edge");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,8.53163006651345e-01-8.53210605848347e-01,"coil +y edge");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,8.54088402519457e-01-8.53210605848347e-01,"coil +z edge");
	pass &= compareResults(MxS->fieldAt(xscan*-1.0)[0]-b0,8.53217812021647e-01-8.53210605848347e-01,"coil -x edge");
	pass &= compareResults(MxS->fieldAt(yscan*-1.0)[0]-b0,8.53163006651345e-01-8.53210605848347e-01,"coil -y edge");
	pass &= compareResults(MxS->fieldAt(zscan*-1.0)[0]-b0,8.54088402519457e-01-8.53210605848347e-01,"coil -z edge");
	
	if(!pass) { vsr::pause(); return pass; }
	printf("Passed bare wires test; building shield...\n");
	
	FieldEstimator2D* fe = new FieldEstimator2D();
	fe->addsource(vec2(-3.92/2,0.61),1.0);
	fe->addsource(vec2(3.92/2,0.61),1.0);
	
	SurfacelCyl* SB = new SurfacelCyl(128);
	SB->retain();
	SB->makeOptCyl(20, 30, .6223, -3.9624/2, 3.9624/2, new PlaneSource(Plane(),10000.0), fe);
	SB->calculateIncident(MxS);
	
	SymmetricSolver SS;
	SS.solve(*SB);
	SS.calculateResult(*SB);
	
	MxS->addsource(SB);
	
	printf("Testing shielded fields...\n");
	b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.60059633351831e+00,"all origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,4.89741926812615e-04,"all +x edge");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-6.08234419476883e-05,"all +y edge");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-1.87974966410431e-04,"all +z edge");
	pass &= compareResults(MxS->fieldAt(xscan*-1.0)[0]-b0,4.89741926807952e-04,"all -x edge");
	pass &= compareResults(MxS->fieldAt(yscan*-1.0)[0]-b0,-6.08234419501308e-05,"all -y edge");
	pass &= compareResults(MxS->fieldAt(zscan*-1.0)[0]-b0,-1.87974966387339e-04,"all -z edge");
	
	vsr::pause();
	MxS->visualize();
	vsr::pause();
	
	MxS->release();
	SB->release();
	
	return pass;
}


bool reference_simpleshield() {
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 3.92);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	
	FieldEstimator2D fe;
	fe.addsource(vec2(-3.92/2,0.55),1.0);
	fe.addsource(vec2(3.92/2,0.55),1.0);
	
	SurfacelCyl* SB = new SurfacelCyl(32);
	SB->retain();
	SB->makeOptCyl(10, 2, .6223, -3.9624/2, 3.9624/2, new PlaneSource(Plane(),10000.0), &fe);
	
	SB->calculateIncident(MxS);
	std::cout << "Incident Response Origin field" << SB->fieldAt(vec3(0,0,0)) << std::endl;
	SB->visualize();
	vsr::pause();
	
	SymmetricSolver SS;
	SS.solve(*SB);
	SS.calculateResult(*SB);
	
	MxS->addsource(SB);
	MxS->visualize();
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);

	bool pass = true;
	mdouble b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.59942943484446e+00,"origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,2.62654994091771e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-3.80699814184204e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-2.60066781668566e-04,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,2.62654994091793e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-3.80699814183982e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,-2.60066781561097e-04,"-z");
	
	vsr::pause();
	SB->release();
	return pass;
}

vec3 f_integ2_test_1(mdouble x, mdouble y, void*) {
	return vec3(x*x + y*y, x + y*y*y, 1 + (x+x*x)*y);
}


bool integrator_tests() {
	Integrator2D I2;
	bool pass = true;
	
	/*
	vec3 v1 = I2.integrate(&f_integ2_test_1, -0.3, 4.6, -6, 7.2);
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);
	*/
	
	return pass;
}

class WavyThing: public DVFunc1<2,mdouble> {
public:
	virtual vec2 operator()(mdouble x) const { return vec2( (x*x-0.5)*3.9624, .6223 + 0.1*(1-x) + 0.04*sin(9*M_PI*x)); }
};

class Ball: public DVFunc1<2,mdouble> {
public:
	Ball(double rr): r(rr) {}
	virtual vec2 operator()(mdouble x) const { return vec2(-r*cos(M_PI*x), r*sin(M_PI*x)); }
	double r;
};


bool csurface_test() {
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 3.92);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	CylSurfaceGeometry SG;
	//SG.zr_profile = new WavyThing;
	Line2D L2D(vec2(-3.9624/2,.6223), vec2(3.9624/2,.6223));
	FieldAdaptiveSurface FAS(L2D);
	SG.zr_profile = &FAS;
	//SG.zr_profile = new Ball(2.2);
	
	SurfaceCurrentRS RS(16,10);
	RS.mySurface = &SG;
	RS.setSurfaceResponse(SurfaceI_Response(10000));
	//RS.setSurfaceResponse(SurfaceI_Response(0));
		
	FAS.optimizeSpacing(fe,0.5);
	RS.calculateIncident(*MxS);
	// Incident Response Origin field< 0.683759 -3.96547e-16 5.57686e-18 >
	// RS.displayContribGrid(vec3(0,0,0),5,5);
	std::cout << "Origin field" << RS.fieldAt(vec3(0,0,0)) << std::endl;
	RS.visualize();
	vsr::pause();
		
	//return true;
	
	SymmetricSolver SS;
	SS.solve(RS);
	SS.calculateResult(RS);
	
	MxS->addsource(&RS);
	
	printf("Testing shielded fields...\n");
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);

	bool pass = true;
	mdouble b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.59942943484446e+00,"origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,2.62654994091771e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-3.80699814184204e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-2.60066781668566e-04,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,2.62654994091793e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-3.80699814183982e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,-2.60066781561097e-04,"-z");
	
	MxS->visualize();
	vsr::pause();
	
	return pass;
}
