#include "tests.hh"
#include "SymmetricSolver.hh"
#include "GenericSolver.hh"
#include "FieldSource.hh"
#include "PlaneSource.hh"
#include "SurfacelCyl.hh"
#include "FieldAnalyzer.hh"
#include "CosThetaBuilder.hh"
#include "Integrator.hh"
#include "SurfaceCurrentSource.hh"
#include "SurfaceCurrentRS.hh"
#include "FieldAdaptiveSurface.hh"
#include "Angles.hh"
#include "UniformField.hh"
#include "SurfaceProfiles.hh"

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

void qsurvey(FieldSource* S, vec3 v1 = vec3(0.5,-0.5,0.5), vec3 v0 = vec3(0,0,0), unsigned int npts = 6) {
	for(float i=0; i<npts; i++) {
		double x = float(i)/(npts-1);
		vec3 l = v1 * x + v0 * (1-x);
		std::cout << "\t" << l << "\t" << S->fieldAt(l) << std::endl;
	}
}

void vis_test_sequence(ReactiveSet* RS, MixedSource* MxS, vec3 vs = vec3(0.5,-0.5,0.5), vec3 vs0 = vec3(0,0,0)) {
	
	std::cout << "B applied: " << std::endl;
	qsurvey(MxS,vs,vs0);
	
	dynamic_cast<MagF_Responder*>(RS)->calculateIncident(*MxS);
	MxS->addsource(dynamic_cast<FieldSource*>(RS));
	MxS->visualize();
	std::cout << "B non-interacting: " << std::endl;
	qsurvey(MxS,vs,vs0);
	vsr::pause();
	
	SymmetricSolver SS;
	SS.solve(*RS);
	SS.calculateResult(*RS);
	MxS->visualize();
	
	std::cout << "B interacting: " << std::endl;
	qsurvey(MxS,vs,vs0);
	
	vsr::pause();
}

void blockade_test() {
	
	MixedSource* MxS = new MixedSource();
	//MxS->loop(-0.25,0.6,3.);
	UniformField* BU = new UniformField(vec3(0,0,1));
	MxS->addsource(BU);
	
	unsigned int nPhi = 8;
	//MagRSCombiner* RSC = new MagRSCombiner(nPhi);
	
	FieldEstimator2Dfrom3D fe(MxS);
	
	Line2D* L2D = new Line2D(vec2(-0.05,1.), vec2(-0.05,0));
	//Dish* D = new Dish(-0.05, -0.25, 1);
	FieldAdaptiveSurface* FAS = new FieldAdaptiveSurface(*L2D);
	FAS->optimizeSpacing(fe, 0.6);
	
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,nPhi,15);
	RS->setSurfaceResponse(SurfaceI_Response(0));
	
	//RSC->addSet(RS);
	
	/*
	Dish* D2 = new Dish(0.5, 0.25, -1);
	CylSurfaceGeometry* SG2 = new CylSurfaceGeometry(D2);
	SurfaceCurrentRS* RS2 = new SurfaceCurrentRS(nPhi,9);
	RS2->mySurface = SG2;
	RS2->setSurfaceResponse(SurfaceI_Response(0));
	RSC.addSet(RS2);
	*/
	
	vis_test_sequence(RS,MxS);
}

void superball_test() {
	
	std::cout << "Expelling external B=(1,1,1) field from a superconducting sphere.\n";
	
	MixedSource* MxS = new MixedSource();
	UniformField* BU = new UniformField(vec3(1,1,1));
	MxS->addsource(BU);
		
	Arc2D* B = new Arc2D(-1.0);
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(B);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,16,15);
	RS->setSurfaceResponse(SurfaceI_Response(0));
	
	vis_test_sequence(RS,MxS);
}


void mirror_test() {
	
	std::cout << "Mirroring cos theta coil current between superconducting plates for improved uniformity.\n";
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.5, 1.0);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_NONE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	
	unsigned int nPhi = 16;
	MagRSCombiner* RSC = new MagRSCombiner(nPhi);

	for(int z=-1; z<=1; z+=2) {
		double w = 0.1;
		RoundedSlab* RSL = new RoundedSlab((0.51+0.5*w)*z, 0.7, w);
		//Line2D* L2D = (z>0)? new Line2D(vec2(0.51, 0.7), vec2(0.51, 0)) : new Line2D(vec2(-0.51, 0), vec2(-0.51, 0.7));
		FieldAdaptiveSurface* FAS = new FieldAdaptiveSurface(*RSL);
		FAS->optimizeSpacing(fe, 0.6);
		
		CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
		SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,nPhi,15);
		RS->setSurfaceResponse(SurfaceI_Response(0));
		RSC->addSet(RS);
	}
	
	vis_test_sequence(RSC,MxS, vec3(0.1,0,0.5), vec3(0.1,0,0));
	
}

void tube_test() {

	std::cout << "Cos theta coil in a ferromagnetic shield.\n";
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 2);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	double t = 0.1;
	double r0 = 0.6;

	RoundedTube* RT = new RoundedTube(vec2(-(1+t),r0+t), vec2(1+t,r0+t), t);
	FieldAdaptiveSurface* FAS = new FieldAdaptiveSurface(*RT);
	FAS->optimizeSpacing(fe, 0.6);
	
	unsigned int nPhi = 16;
	
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,nPhi, 17);
	RS->setSurfaceResponse(SurfaceI_Response(10000));
	vis_test_sequence(RS, MxS, vec3(0.1,0,0.5), vec3(0.1,0,0));
}


void flux_trap_test() {

	std::cout << "Calculating trapped-flux state of superconducting ring..." << std::endl;
	
	RoundedTube* RT = new RoundedTube(vec2(0,0.25), vec2(0,0.8), 0.1);
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(RT);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, 10, 18);
	RS->setSurfaceResponse(SurfaceI_Response(0));
	RS->visualize();
	
	SymmetricSolver SS;
	SS.solve(*RS);
	
	std::cout << "Setting initial current loop distribution..." << std::endl;
	RS->setZeroState();
	//RS->set_current_loop(0,1.0,false);
	RS->set_current_loop(9,1.0);
	RS->setFinalState(RS->finalState.normalized()*(0.05*RS->nDF()));
	RS->visualize();
	vsr::pause();
		
	for(unsigned int i=0; i<100; i++) {
		
		std::cout << "Iterating solution to stable trapped-flux state..." << std::endl;
		RS->incidentState = RS->finalState;
		//SS.calculateResult(*RS);	// alternate, gets there faster, but with huge normalization isse
		SS.selfInteract(*RS);
		
		std::cout << RS->finalState.mag() << std::endl;
		RS->setFinalState(RS->finalState.normalized()*(0.05*RS->nDF()));
		std::cout << "Initial/final state difference magnitude: " << (RS->finalState - RS->incidentState).mag() << std::endl;
		
		qsurvey(RS, vec3(0,.75,0), vec3(0,0,0), 6);
		
		RS->visualize();
		vsr::pause();
	}

}

bool csurface_test() {

	// mimic "simple shield" test with continuous surface

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
	
	Line2D L2D(vec2(-zh,r0), vec2(zh,r0));
	FieldAdaptiveSurface FAS(L2D);
	FAS.optimizeSpacing(sfe,0.3, false);
	FAS.symmetry_test();
	CylSurfaceGeometry SG(&FAS);
	SurfaceCurrentRS RS(&SG,nPhi,12);
	RS.setSurfaceResponse(SurfaceI_Response(10000));
	RSC.addSet(&RS);
	
	
	RSC.calculateIncident(*MxS);

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
	SurfaceCurrentRS RS(&SG,128,50);
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
