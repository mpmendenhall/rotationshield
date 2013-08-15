/// file "main.cpp" \brief Menu system for launching various shield simulations

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#include "CosThetaBuilder.hh"
#include "SymmetricSolver.hh"
#include "PlaneSource.hh"
#include "ShieldBuilder.hh"
#include "analysis.hh"
#include "PathUtils.hh"
#include "ControlMenu.hh"
#include "tests.hh"
#include "Studies.hh"


void mi_sampleShieldVis(std::deque<std::string>&, std::stack<std::string>&) {
	// set up output paths
	std::string basepath = "../";
	std::string projectpath = basepath+"/CosThetaCoil/";
	makePath(projectpath);
	std::string fieldspath = projectpath+"/Fields/";
	makePath(fieldspath);
	std::string gfpath = projectpath+"/GFs/";
	makePath(gfpath);
	
	// set up output files for field scans
	std::ofstream fieldsout;
	std::ofstream statsout;
	fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
	statsout.open((fieldspath+"/Fieldstats.txt").c_str());
	
	// set up cos theta coil
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.648;	// coil radius
	mdouble clen = 4.292;		// coil length
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;	// coil distortion parameter
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner);
	b.regularCoil(*ms,&pert);
	ms->visualize();
	//vsr::pause();
	
	// optional: construct shield
	if(1) {
		mdouble shieldDistance = 0.069;				// distance of shield from wires, 6.9cm
		mdouble sradius = cradius + shieldDistance;	// shield radius
		mdouble slen = clen + 0.40;					// shield length, coil + 40cm
		FieldEstimator2D* fe = new FieldEstimator2D();
		fe->addsource(vec2(-clen/2.0,cradius),1.0);
		fe->addsource(vec2(clen/2.0,cradius),1.0);
		ShieldBuilder* G = new ShieldBuilder(128);	// 128 segments around shield builder
		G->makeOptCyl(10, 20, sradius, -slen/2, slen/2, new PlaneSource(Plane(),10000), fe); // optimized grid with 10+20 divisions along z
		// solve for shield response function
		SymmetricSolver* sp = new SymmetricSolver(G);
		sp->solve();
		// calculate response to coil fields
		sp->calculateIncident(ms);
		sp->calculateResult();
		ms->addsource(sp);
		ms->visualize();
	}
	
	// example of probing field value
	vec3 b0 = ms->fieldAt(vec3(0,0,0));
	printf("Field at center: ");
	b0.display();
	
	// survey data points on grid
	FieldAnalyzer FA = FieldAnalyzer(ms);
	//        lower left corner,     upper right corner, nx,ny,nz, output files
	FA.survey(vec3(0.0,0.0,0.0),vec3(0.20,0.20,0.5),10,10,10,statsout,fieldsout);
	
	
	// cleanup
	ms->release();
	fieldsout.close();
	statsout.close();
	vsr::pause();
}

void mi_runtests(std::deque<std::string>&, std::stack<std::string>&) {
	printf("Simple shield self-test...\n");
	reference_simpleshield();
	printf("\n\nReference self-test...\n");
	reference_sanity_check();
}

void mi_ShortCoilOptRad(std::deque<std::string>&, std::stack<std::string>& stack) {
	
	float rd = streamInteractor::popFloat(stack);
	printf("Generating short coil simulation for r=%g m\n",rd);
	
	// set up output paths and files
	std::string projectpath = getEnvSafe("ROTSHIELD_OUT")+"/HalfCoil_"+dtos(rd);
	std::string fieldspath = projectpath+"/Fields/";
	makePath(fieldspath);
	std::ofstream fieldsout;
	std::ofstream statsout;
	fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
	statsout.open((fieldspath+"/Fieldstats.txt").c_str());
	
	// run
	shortCoil(fieldsout,statsout,rd);

	// cleanup
	fieldsout.close();
	statsout.close();
}


void menuSystem(std::deque<std::string> args=std::deque<std::string>()) {
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester peek("Show stack",&menutils_PrintStack);
	
	inputRequester sampleShieldVis("Example shield with visualization",&mi_sampleShieldVis);
	inputRequester selfTests("Self test to reproduce known results",&mi_runtests);
	inputRequester shortCoil("'Vertical' Short coil",&mi_ShortCoilOptRad);
	shortCoil.addArg("Radius","1.0");
		
	OptionsMenu OM("Rotation Shield Main Menu");
	//OM.addChoice(&sampleShieldVis,"exm");
	//OM.addChoice(&selfTests,"test");
	OM.addChoice(&shortCoil,"short");
	OM.addChoice(&exitMenu,"x");
	OM.addChoice(&peek,"peek",SELECTOR_HIDDEN);
	std::stack<std::string> stack;
	OM.doIt(args,stack);
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
}

void* menuThread(void* args) {
	std::deque<std::string>& inArgs = *(std::deque<std::string>*)args;
	menuSystem(inArgs);
	return NULL;
}

int main(int argc, char *argv[]) {
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
#ifdef WITH_OPENGL
	vsr::initWindow();
	pthread_t thread;
	pthread_create( &thread, NULL, &menuThread, &args );
	vsr::doGlutLoop();
#else
	menuSystem(args);
#endif
	return 0;
}
