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
#include "Angles.hh"


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
	SB->calculateIncident(*MxS);
	
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
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
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
	
	// non-interacting field reaction
	SB->calculateIncident(*MxS);
	bool pass = true;
	mdouble b0 = SB->fieldAt(origin)[0];
	pass &= compareResults(b0,6.83758710057794e-01,"origin - noninteracting");
	pass &= compareResults(SB->fieldAt(xscan)[0]-b0,1.11572338749721e-03,"+x");
	pass &= compareResults(SB->fieldAt(yscan)[0]-b0,-9.86029434321134e-05,"+y");
	pass &= compareResults(SB->fieldAt(zscan)[0]-b0,-1.33717928505817e-03,"+z");
	pass &= compareResults(SB->fieldAt(-xscan)[0]-b0,1.11572338749666e-03,"-x");
	pass &= compareResults(SB->fieldAt(-yscan)[0]-b0,-9.86029434315583e-05,"-y");
	pass &= compareResults(SB->fieldAt(-zscan)[0]-b0,-1.33717928543731e-03,"-z");
	SB->visualize();
	vsr::pause();
	
	SymmetricSolver SS;
	SS.solve(*SB);
	SS.calculateResult(*SB);
	
	MxS->addsource(SB);
	MxS->visualize();

	b0 = MxS->fieldAt(origin)[0];
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

mvec f_integ2_test_1(vec2 xy, void*) {
	double x = xy[0];
	double y = xy[1];
	return mvec(vec3(x*x + y*y, x + y*y*y, 1 + (x+x*x)*y));
}


bool integrator_tests() {
	Integrator2D I2;
	I2.setMethod(INTEG_GSL_QAG);
	
	bool pass = true;
		
	mvec v1 = I2.integrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2));
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);

	v1 = I2.polarIntegrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2), vec2(-4,6));
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);
	
	v1 = I2.polarIntegrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2), vec2(2,1));
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);
	
	return pass;
}

void angular_intervals_test() {
	
	Angular_Interval_Set AIS;
	
	AIS.add_interval(angular_interval(0.1*M_PI, 0.2*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.3*M_PI, 0.4*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.2*M_PI, 0.25*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.25*M_PI, 0.3*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.5*M_PI, 0.6*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.42*M_PI, 0.48*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.add_interval(angular_interval(0.39*M_PI, 0.51*M_PI));
	std::cout << AIS << std::endl;
	
	std::cout << std::endl;
	
	AIS.subtract_interval(angular_interval(0.3*M_PI, 0.5*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.subtract_interval(angular_interval(0.2*M_PI, 0.3*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.subtract_interval(angular_interval(0.1*M_PI, 0.21*M_PI));
	std::cout << AIS << std::endl;
	
	AIS.subtract_interval(angular_interval(-0.1*M_PI, M_PI));
	std::cout << AIS << std::endl;
	
	std::cout << std::endl;
	
	AIS.add_interval(angular_interval(-0.5*M_PI, 0.5*M_PI));
	std::cout << AIS << std::endl;
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

	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 3.92);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	FieldEstimator2D sfe;
	sfe.addsource(vec2(-3.92/2,0.61),1.0);
	sfe.addsource(vec2(3.92/2,0.61),1.0);
	
	unsigned int nPhi = 32;
	
	MagRSCombiner RSC(nPhi);
	
	double zh = 3.9624/2;
	double r0 = .6223;
	
	// main shield
	Line2D L2D(vec2(-zh,r0), vec2(zh,r0));
	FieldAdaptiveSurface FAS(L2D);
	FAS.optimizeSpacing(sfe,0.3, false);
	FAS.symmetry_test();
	CylSurfaceGeometry SG(&FAS);
	SurfaceCurrentRS RS(nPhi,12);
	RS.mySurface = &SG;
	RS.setSurfaceResponse(SurfaceI_Response(10000));
	RSC.addSet(&RS);
	
	/*
	// rear SC endcap 
	Line2D L_Endcap(vec2(-zh,0), vec2(-zh,r0));
	FieldAdaptiveSurface FAS_EC(L_Endcap);
	FAS_EC.optimizeSpacing(fe,0.5);
	CylSurfaceGeometry SG_EC(&FAS_EC);
	SurfaceCurrentRS RS_EC(nPhi,12,0);
	RS_EC.mySurface = &SG_EC;
	RS_EC.setSurfaceResponse(SurfaceI_Response(0));
	RSC.addSet(&RS_EC);
	*/
	
	RSC.calculateIncident(*MxS);
	std::cout << "Net shield current: " << RSC.net_current() << std::endl;
	
	bool pass = true;
	mdouble b0 = RSC.fieldAt(origin)[0];
	pass &= compareResults(b0,6.83758710057794e-01,"origin - noninteracting");
	pass &= compareResults(RSC.fieldAt(xscan)[0]-b0,1.11572338749721e-03,"+x");
	pass &= compareResults(RSC.fieldAt(yscan)[0]-b0,-9.86029434321134e-05,"+y");
	pass &= compareResults(RSC.fieldAt(zscan)[0]-b0,-1.33717928505817e-03,"+z");
	pass &= compareResults(RSC.fieldAt(-xscan)[0]-b0,1.11572338749666e-03,"-x");
	pass &= compareResults(RSC.fieldAt(-yscan)[0]-b0,-9.86029434315583e-05,"-y");
	pass &= compareResults(RSC.fieldAt(-zscan)[0]-b0,-1.33717928543731e-03,"-z");

	RSC.visualize();
	vsr::pause();
		
	//return true;
	
	SymmetricSolver SS;
	SS.solve(RSC);
	SS.calculateResult(RSC);
	
	MxS->addsource(&RSC);
	
	std::cout << "Net shield current: " << RSC.net_current() << std::endl;
	printf("Testing shielded fields...\n");

	b0 = MxS->fieldAt(origin)[0];
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

bool csurface_test_B() {
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	
	CosThetaBuilder CB(17, 0.61, 3.92);
	CB.myCap[0] = CB.myCap[1] = CosThetaBuilder::CAP_LINE;
	CB.buildCoil(*MxS);
	MxS->visualize();
	
	FieldEstimator2D* fe = new FieldEstimator2D();
	fe->addsource(vec2(-3.92/2,0.61),1.0);
	fe->addsource(vec2(3.92/2,0.61),1.0);

	// main shield
	double zh = 3.9624/2;
	double r0 = .6223;
	MagRSCombiner RSC(128);
	Line2D L2D(vec2(-zh,r0), vec2(zh,r0));
	FieldAdaptiveSurface FAS(L2D);
	FAS.optimizeSpacing(*fe,0.2);
	CylSurfaceGeometry SG(&FAS);
	SurfaceCurrentRS RS(128,50);
	RS.mySurface = &SG;
	RS.setSurfaceResponse(SurfaceI_Response(10000));
	RSC.addSet(&RS);

	RSC.calculateIncident(*MxS);
	std::cout << "Origin field" << RSC.fieldAt(vec3(0,0,0)) << std::endl;
	RSC.visualize();
	vsr::pause();
			
	SymmetricSolver SS;
	SS.solve(RSC);
	SS.calculateResult(RSC);
	
	MxS->addsource(&RSC);

	printf("Testing shielded fields...\n");
	bool pass = true;
	double b0 = MxS->fieldAt(origin)[0];
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

	return pass;
}
