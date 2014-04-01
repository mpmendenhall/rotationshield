/* 
 * Tests.cpp, part of the RotationShield program
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
#include "PathUtils.hh"
#include "GenericSolver.hh"
#include "DipoleSource.hh"
#include "HoleDipolePerturbation.hh"
#include "LineSource.hh"
#include "MultiQuilibrator.hh"

#include <ctime>


bool compareResults(double a, double b, const char* label) {
	bool pass = true;
	if(label) printf("%s:\n",label);
	printf("\tTest Result:\t%.14e\n\t  Reference:\t%.14e\n\t\t   error:\t%.3f%%",double(a),double(b),double(100*(a-b)/a));
	if(fabs(100*(a-b)/a) > 5.0 ) { printf(" <<-------- ******* TEST FAILED ******"); pass = false; }
	printf("\n");
	fflush(stdout);
	return pass;
}

bool reference_simpleshield() {

	// Simplified coil in crude ferromagnetic shield
	// test guards against changes breaking basic line source field measurements or SymmetricSolver algoithm
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
	MagExtField* MxS = new MagExtField();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 3.92);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	
	// Bare wires field
	bool pass = true;
	double b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,9.48809810037848e-01,"origin - bare wires");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,1.50039293277437e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-2.76705054540693e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,1.01000713313459e-03,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,1.50039293277426e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-2.76705054540805e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,1.01000713313459e-03,"-z");

	if(!pass) { return pass; }
		  
	FieldEstimator2D fe;
	fe.addsource(vec2(-3.92/2,0.55),1.0);
	fe.addsource(vec2(3.92/2,0.55),1.0);
	
	SurfacelCyl* SB = new SurfacelCyl(32);
	SB->retain();
	SB->makeOptCyl(10, 2, .6223, -3.9624/2, 3.9624/2, new PlaneSource(Plane(),10000.0), &fe);
	
	// non-interacting field reaction
	SB->incidentState = SB->getFullReactionTo(MxS);
	b0 = SB->fieldAt(origin)[0];
	pass &= compareResults(b0,6.83758710057794e-01,"origin - noninteracting shield only");
	pass &= compareResults(SB->fieldAt(xscan)[0]-b0,1.11572338749721e-03,"+x");
	pass &= compareResults(SB->fieldAt(yscan)[0]-b0,-9.86029434321134e-05,"+y");
	pass &= compareResults(SB->fieldAt(zscan)[0]-b0,-1.33717928505817e-03,"+z");
	pass &= compareResults(SB->fieldAt(-xscan)[0]-b0,1.11572338749666e-03,"-x");
	pass &= compareResults(SB->fieldAt(-yscan)[0]-b0,-9.86029434315583e-05,"-y");
	pass &= compareResults(SB->fieldAt(-zscan)[0]-b0,-1.33717928543731e-03,"-z");
	SB->visualize();
	vsr::pause();
	
	if(!pass) { return pass; }
	
	SymmetricSolver SS;
	SS.solve(*SB);
	SS.calculateResult(*SB);
	
	MxS->addsource(SB);
	MxS->visualize();
	
	// interacting field reaction
	b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.59941232856854e+00,"origin, wires + interacting shield");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,2.62856023674884e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-3.80588069445853e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-2.67641333622226e-04,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,2.62856023674840e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-3.80588069445409e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,-2.67641333894009e-04,"-z");
	
	vsr::pause();
	SB->release();
	return pass;
}

bool reference_simpleshield_cached() {

	// Test solution cache IO
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
	MagExtField* MxS = new MagExtField();
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
	SB->incidentState = SB->getFullReactionTo(MxS);
	
	SymmetricSolver* SS = SymmetricSolver::cachedSolve(*SB,getEnvSafe("ROTSHIELD_OUT",".")+"/ref_simpleshield_sol.dat");
	SS->calculateResult(*SB);
	
	MxS->addsource(SB);
	MxS->visualize();
	
	// interacting field reaction
	bool pass = true;
	double b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.59941232856854e+00,"origin, wires + interacting shield");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,2.62856023674884e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-3.80588069445853e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-2.67641333622226e-04,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,2.62856023674840e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-3.80588069445409e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,-2.67641333894009e-04,"-z");
	
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
	IntegratorND IN;
	
	bool pass = true;
		
	mvec v1 = I2.integrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2));
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);

	v1 = I2.polarIntegrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2), vec2(-4,6));
	pass &= compareResults(v1[0], 1390.8356);
	pass &= compareResults(v1[1], 1843.50936);
	pass &= compareResults(v1[2], 405.15552);
	
	//v1 = I2.polarIntegrate2D(&f_integ2_test_1, vec2(-0.3,-6), vec2(4.6, 7.2), vec2(2,1));
	v1 = IN.integratePolar(&f_integ2_test_1, 3, vec2(2,1), vec2(-0.3,-6), vec2(4.6, 7.2));
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


void vis_test_sequence(ReactiveSet* RS, MagExtField* MxS, vec3 vs = vec3(0.5,-0.5,0.5), vec3 vs0 = vec3(0,0,0)) {
	
	std::cout << "B applied: " << std::endl;
	qsurvey(MxS,vs,vs0);
	
	RS->incidentState = RS->getFullReactionTo(MxS);
	RS->setFinalState(RS->incidentState);
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


void superball_test(unsigned int ngrid) {
	
	std::cout << "Expelling external B=(1,1,1) field from a superconducting sphere.\n";
	
	MagExtField* MxS = new MagExtField();
	UniformField* BU = new UniformField(vec3(1,1,1));
	MxS->addsource(BU);
		
	Arc2D* B = new Arc2D(-1.0);
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(B);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,ngrid,ngrid);
	RS->setSurfaceResponse(SurfaceI_Response(0));
	
	vis_test_sequence(RS,MxS);
}


void mirror_test() {
	
	std::cout << "Mirroring cos theta coil current between superconducting plates for improved uniformity.\n";
	
	MagExtField* MxS = new MagExtField();
	CosThetaBuilder b = CosThetaBuilder(5, 0.5, 1.0);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_NONE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	
	unsigned int nPhi = 16;
	MagRSCombiner* RSC = new MagRSCombiner(nPhi);

	for(int z=-1; z<=1; z+=2) {
		double w = 0.1;
		RoundedSlab* RSL = new RoundedSlab((0.51+0.5*w)*z, 0.7, w);
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
	
	MagExtField* MxS = new MagExtField();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 2);
	b.myCap[0] = b.myCap[1] = CosThetaBuilder::CAP_LINE;
	b.buildCoil(*MxS);
	FieldEstimator2Dfrom3D fe(MxS);
	
	double t = 0.1;
	double r0 = 0.6;

	RoundedTube* RT = new RoundedTube(vec2(-(1+t),r0+t), vec2(1+t,r0+t), t);
	FieldAdaptiveSurface* FAS = new FieldAdaptiveSurface(*RT);
	FAS->optimizeSpacing(fe, 0.6);
	
	unsigned int nPhi = 32;
	
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG,nPhi, 25);
	RS->setSurfaceResponse(SurfaceI_Response(10000));
	vis_test_sequence(RS, MxS, vec3(0.1,0,0.5), vec3(0.1,0,0));
}

void flux_trap_test() {

	std::cout << "Calculating trapped-flux state of superconducting ring..." << std::endl;
	
	RoundedTube* RT = new RoundedTube(vec2(-.2, 0.25), vec2(0.2, 0.8), 0.1);
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(RT);
	SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, 10, 18);
	RS->setSurfaceResponse(SurfaceI_Response(0));
	RS->visualize();
	
	SymmetricSolver SS;
	SS.solve(*RS);
	
	std::cout << "Setting initial current loop distribution..." << std::endl;
	RS->setZeroState();
	RS->set_current_loop(9,1.0);
	RS->setFinalState(RS->finalState.normalized()*(0.05*RS->nDF()));
	std::cout << "Total state vector magnitude: " << RS->finalState.mag() << std::endl;
	RS->visualize();
	vsr::pause();
	
	std::cout << "Selecting singular component from SVD, which indicates ``self-interacting'' mode without external field:" << std::endl;
	RS->incidentState = RS->finalState;
	SS.set_singular_epsilon(-0.01);	// pick out only singular component of solution
	SS.calculateResult(*RS);		// apply singular state
	RS->visualize();
	std::cout << "Total state vector magnitude: " << RS->finalState.mag() << std::endl;
	qsurvey(RS, vec3(0,.75,0), vec3(0,0,0), 6);
	vsr::pause();
	
	for(unsigned int i=0; i<5; i++) {
		
		std::cout << "Checking stability of trapped-flux state by repeatedly applying interaction operator..." << std::endl;
		RS->incidentState = RS->finalState;
		SS.selfInteract(*RS);	// state should be stable under self-interaction
		std::cout << "Total state vector magnitude: " << RS->finalState.mag() << std::endl;
		std::cout << "Initial/final state difference magnitude: " << (RS->finalState - RS->incidentState).mag() << std::endl;
		
		qsurvey(RS, vec3(0,.75,0), vec3(0,0,0), 6);
		
		RS->visualize();
		vsr::pause();
	}
}

void hole_perturbation_test() {

	Arc2D* B = new Arc2D(-0.7);
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(B);
	
	MagRSCombiner Holes(1);
	for(unsigned int i=0; i<10; i++)
		Holes.addSet(new HoleDipolePerturbation(*SG, vec2(randunif(0,1),randunif(0,1)), 0.1));
	
	Holes.visualize();
	
		
	GenericSolver GS;
	GS.solve(Holes);
	vsr::pause();
}

void two_torus_equilibration() {

	// set up interlinking torus geometry
	Arc2D* B = new Arc2D(0.2, -M_PI, M_PI, vec2(0,0.6));
	CylSurfaceGeometry* SG = new CylSurfaceGeometry(B);
	TransformedGeometry* TG1 = new TransformedGeometry(SG);
	TransformedGeometry* TG2 = new TransformedGeometry(SG);
	TG1->transl_v = vec3(0,0.3,0);
	TG2->transl_v = vec3(0,-0.3,0);
	TG1->transf_M = Matrix<3,3,double>::rotation(0,2,0.14*M_PI);
	TG2->transf_M = Matrix<3,3,double>::rotation(0,2,-0.14*M_PI);
	
	// create ReactiveSets with different number of DF
	SurfaceCurrentRS* RS1 = new SurfaceCurrentRS(TG1, 16, 8);
	RS1->setSurfaceResponse(SurfaceI_Response(0));
	SurfaceCurrentRS* RS2 = new SurfaceCurrentRS(TG2, 8, 7);
	RS2->setSurfaceResponse(SurfaceI_Response(0));
	// combine them into one unit
	MagRSCombiner* MRC = new MagRSCombiner(1);
	MRC->addSet(RS1);
	MRC->addSet(RS2);
	MRC->visualize();
	
	// external incident field source
	MagExtField* MxI = new MagExtField();
	UniformField* BU = new UniformField(vec3(.1,.1,.3));
	MxI->addsource(BU);
	// combined external and surface fields
	MixedSource* Ftot = new MixedSource();
	Ftot->addsource(MxI);
	Ftot->addsource(MRC);
	
	// set incident state
	MRC->incidentState = MRC->getFullReactionTo(MxI); // TODO make this propagate down
	RS1->incidentState = RS1->getFullReactionTo(MxI);
	RS2->incidentState = RS2->getFullReactionTo(MxI);
	std::cout << "Initial field RMS inside 1: " << MxI->field_RMS_near(*RS1->mySurface,-0.1) << "\n";
	std::cout << "Initial field RMS inside 2: " << MxI->field_RMS_near(*RS2->mySurface,-0.1) << "\n";
	
	// start timer
	time_t start_time, end_time;
	start_time = std::time(NULL);
	
	// solve each individually
	SymmetricSolver SS1;
	SS1.solve(*RS1);
	SS1.set_singular_epsilon(1e-4);
	//
	SymmetricSolver SS2;
	SS2.solve(*RS2);
	SS2.set_singular_epsilon(1e-4);
	
	SS1.calculateResult(*RS1);
	SS2.calculateResult(*RS2);
	MRC->load_component_DF();
	//
	MRC->visualize();
	//std::cout << "Non-interacting field RMS inside 1: " << Ftot->field_RMS_near(*RS1->mySurface,-0.1) << "\n";
	//std::cout << "Non-interacting field RMS inside 2: " << Ftot->field_RMS_near(*RS2->mySurface,-0.1) << "\n";
	
	// iterate to stability
	MultiQuilibrator MQ;
	MQ.addSet(RS1,&SS1);
	MQ.addSet(RS2,&SS2);
	while(MQ.step() > 1e-4)
		MRC->visualize();
	
	// end timer
	end_time = std::time(NULL);
	
	std::cout << "Interacting field RMS inside 1: " << Ftot->field_RMS_near(*RS1->mySurface,-0.1) << "\n";
	std::cout << "Interacting field RMS inside 2: " << Ftot->field_RMS_near(*RS2->mySurface,-0.1) << "\n";
	
	std::cout << "\nIterative method calculation complete in " << difftime(end_time,start_time) << " seconds.\n";
	std::cout << "Calculating brute-force solution...\n\n";
	
	MRC->load_component_DF();
	mvec final_iterative = MRC->finalState;
	
	////////////////
	// Direct brute-force calculation
	
	// start timer
	start_time = std::time(NULL);
	
	// solve
	GenericSolver GS;
	GS.solve(*MRC);
	GS.set_singular_epsilon(2e-3);
	GS.calculateResult(*MRC);
	
	// end timer
	end_time = std::time(NULL);
	
	std::cout << "Interacting field RMS inside 1: " << Ftot->field_RMS_near(*RS1->mySurface,-0.1) << "\n";
	std::cout << "Interacting field RMS inside 2: " << Ftot->field_RMS_near(*RS2->mySurface,-0.1) << "\n";
	std::cout << "\nBrute-force calculation complete in " << difftime(end_time,start_time) << " seconds.\n";
	std::cout << "Final state difference = " << (final_iterative - MRC->finalState).mag()/MRC->finalState.mag() << "%\n";
	
	vsr::pause();
	MRC->visualize();
}

