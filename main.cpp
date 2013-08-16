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

/// run self-tests
void mi_runtests(std::deque<std::string>&, std::stack<std::string>&) {
	printf("Simple shield self-test...\n");
	reference_simpleshield();
	printf("\n\nReference self-test...\n");
	reference_sanity_check();
}

/// Short ("vertical") coil with re-oriented sample cell
/// N=30, length=2.m, radius variable
/// shield radius nominal 10cm outside coil; shield 40cm longer than coil
/// cell dimensions y=40cm, x=7.5cm, z=10cm, inner edge 5cm from center in x (field) direction
/// Want:
///   B0 = 30mG; <dBi/di> < 10^-7 G/cm
///   sqrt(<(dBx/dz)^2>) (z=long direction, x=along polarization)
///   cells: | 7.5 | 10 | 7.5 |    --> x
void mi_ShortCoilOptRad(std::deque<std::string>&, std::stack<std::string>& stack) {
	float rd = streamInteractor::popFloat(stack);
	printf("Generating short coil simulation for r=%g m\n",rd);
	std::string projectpath = getEnvSafe("ROTSHIELD_OUT")+"/ShortCoil_"+dtos(rd);
	
	nEDM_Geom GM(projectpath);
	GM.coil = new CosThetaBuilder(15,rd,2.0);
	GM.shield = new coilShield();
	GM.shield->length = GM.coil->length+0.4;
	GM.shield->radius = GM.coil->radius+0.1;
	GM.construct();
	
	GM.cell = new fieldCell(vec3(0.05,-0.20,-0.05),vec3(0.125,0.20,0.05),9,31,9);
	GM.takeSample();
	
	printf("Data collection complete.\n");
}


/// Bare coil, variable N/r/l geometry
void mi_BareCoil(std::deque<std::string>&, std::stack<std::string>& stack) {
	float rd = streamInteractor::popFloat(stack);
	float ln = streamInteractor::popFloat(stack);
	int hn = streamInteractor::popInt(stack);
	
	printf("Generating coil simulation for n=%i, l=%g, r=%g m\n",hn,ln,rd);
	std::string projectpath = getEnvSafe("ROTSHIELD_OUT")+"/BareCoil_"+itos(hn)+"_"+dtos(ln)+"_"+dtos(rd);
	
	nEDM_Geom GM(projectpath);
	GM.coil = new CosThetaBuilder(hn,rd,ln);
	GM.construct();
	GM.cell = new fieldCell(vec3(-0.2,-0.2,-0.2),vec3(0.2,0.2,0.2),41,41,41);
	GM.takeSample();
	
	printf("Data collection complete.\n");
}


void menuSystem(std::deque<std::string> args=std::deque<std::string>()) {
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester peek("Show stack",&menutils_PrintStack);
	
	inputRequester selfTests("Self test to reproduce known results",&mi_runtests);
	
	inputRequester shortCoil("'Vertical' Short coil",&mi_ShortCoilOptRad);
	shortCoil.addArg("Radius","1.0");
	
	inputRequester bareCoil("Generic bare coil field probe",&mi_BareCoil);
	bareCoil.addArg("half n","15");
	bareCoil.addArg("Length","4.292");
	bareCoil.addArg("Radius","0.648");
	
	OptionsMenu OM("Rotation Shield Main Menu");
	//OM.addChoice(&selfTests,"test");
	OM.addChoice(&shortCoil,"short");
	OM.addChoice(&bareCoil,"bare");
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
