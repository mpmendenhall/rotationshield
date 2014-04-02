/// file "main.cpp" \brief Menu system for launching various shield simulations
/* Solver for magnetic fields in rotationally symmetric configurations of linear media
 *
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * this project re-distributes "cubature," (c) 2005-2013 Steven G. Johnson,
 * distributed under GPL v2 or later
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

#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#include "PathUtils.hh"
#include "Typedefs.hh"
#include "ControlMenu.hh"
#include "tests.hh"
#include "Studies.hh"
#include "ncube.hh"

template<unsigned int N>
void hypercube_rotator() {
	NRotator<N> NR;
	NCube<N> NC;
	
	double ftime = 15.e-3; // min frame time in s
	srand(time(NULL));
	std::vector<double> vrot;
	for(unsigned int i=1; i<N; i++)
		for(unsigned int j=0; j<i; j++)
			vrot.push_back(randunif(-2./sqrt(N),2./sqrt(N)));
	
	vsr::set_pause();
	if(vsr::get_pause())
		printf("Press [ENTER] in visualization window to continue...\n");
		
	clock_t prevTime = clock();
	while(vsr::get_pause()) {
		clock_t timeNow = clock();
		double dtime = double(timeNow-prevTime)/double(CLOCKS_PER_SEC);
		prevTime = timeNow;
		unsigned int n=0;
		for(unsigned int i=1; i<N; i++)
			for(unsigned int j=0; j<i; j++)
				NR.rotate(i, j, vrot[n++]*dtime);
		NC.visualize(&NR);
		if(dtime>=0 && dtime<ftime)
			usleep((ftime-dtime)*1.e6);
	}
}

/// run self-tests
void mi_runtests(StreamInteractor*) {
	hole_perturbation_test();
	return;
	printf("Simple shield self-test...\n");
	reference_simpleshield();
}

void mi_demo_sphere(StreamInteractor* S) {
	int nGrid = S->popInt();
	if(nGrid < 4 || nGrid > 256) { printf("Are you kidding?\n"); return; }
	superball_test(nGrid);
}
void mi_demo_mirror(StreamInteractor*) { mirror_test(); }
void mi_demo_tube(StreamInteractor*) { tube_test(); }
void mi_demo_ring(StreamInteractor*) { flux_trap_test(); }
void mi_demo_tori(StreamInteractor*) { two_torus_equilibration(); }

/// hypercube visualization test
void mi_ncube(StreamInteractor* S) {
	int n = S->popInt();
	if(n==3) hypercube_rotator<3>();
	else if(n==4) hypercube_rotator<4>();
	else if(n==5) hypercube_rotator<5>();
	else printf("Only 3 to 5 dimensions allowed.\n");
}

/// set clear color for visualization
void mi_clearcolor(StreamInteractor* S) {
	float b = S->popFloat();
	float g = S->popFloat();
	float r = S->popFloat();
	vsr::setClearColor(r,g,b);
	vsr::startRecording();
	vsr::clearWindow();
	vsr::stopRecording();
}

void menuSystem(std::deque<std::string> args = std::deque<std::string>()) {
	
	InputRequester ncube("Hyercube visualization test", &mi_ncube);
	ncube.addArg("n dim","3");
	InputRequester setClearColor("Set visualization background color", &mi_clearcolor);
	setClearColor.addArg("r","0.0");
	setClearColor.addArg("g","0.0");
	setClearColor.addArg("b","0.0");
	
	InputRequester selfTests("Self test to reproduce known results", &mi_runtests);
	
	OptionsMenu selectDemo("Demonstration");
	InputRequester demo_sphere("superconducting sphere",&mi_demo_sphere);
	demo_sphere.addArg("gridding","16");
	InputRequester demo_tube("cos-theta coil in ferromagnetic tube",&mi_demo_tube);
	InputRequester demo_ring("trapped-flux state of superconducting ring",&mi_demo_ring);
	InputRequester demo_mirror("superconductor-mirrored cos theta coil",&mi_demo_mirror);
	InputRequester demo_tori("interlinked superconducting rings",&mi_demo_tori);
	
	selectDemo.addChoice(&demo_sphere,"sc");
	selectDemo.addChoice(&demo_mirror,"mr");
	selectDemo.addChoice(&demo_tube,"tb");
#ifdef WITH_LAPACKE
	selectDemo.addChoice(&demo_ring,"ft");
	selectDemo.addChoice(&demo_tori,"ir");
#endif
	selectDemo.addChoice(&InputRequester::exitMenu,"x");

	
	OptionsMenu OM("Rotation Shield Main Menu");
#ifdef WITH_OPENGL
	OM.addChoice(&ncube,"ncube",SELECTOR_HIDDEN);
	OM.addChoice(&setClearColor,"cc",SELECTOR_HIDDEN);
#endif
	OM.addChoice(&selectDemo,"demo");
	
	SystemConfiguration SC;
	OM.addChoice(&SC.outDir,"dir");
	OM.addChoice(&SC.OMfieldsrc,"field");
	OM.addChoice(&SC.OMsurfaces,"bound");
	OM.addChoice(&SC.OMcell,"cell");
	
	OM.addChoice(&SC.doSolve,"solve");
	OM.addChoice(&SC.doApply,"apply");
	OM.addChoice(&SC.equilibratePtb,"ptb",SELECTOR_HIDDEN);
	OM.addChoice(&SC.zeroResponse,"zero");
	OM.addChoice(&SC.doMeas,"meas");
	OM.addChoice(&SC.qSurvey,"qf");
#ifdef WITH_LAPACKE
	OM.addChoice(&SC.addSingular,"svd");
	OM.addChoice(&SC.setSingularEpsilon,"ep");
#endif
	
	OM.addChoice(&InputRequester::exitMenu,"x");
	OM.addChoice(&selfTests,"test",SELECTOR_HIDDEN);
	
	std::stack<std::string> stack;
	OM.mydeque = &args;
	OM.mystack = &stack;
	OM.doIt();
	
	printf("\n\n\n>>>>> Goodbye. <<<<<\n\n\n");
}

void* menuThread(void* args) {
	std::deque<std::string>& inArgs = *(std::deque<std::string>*)args;
	menuSystem(inArgs);
	vsr::set_kill();
	return NULL;
}

int main(int argc, char *argv[]) {
	std::deque<std::string> args;
	for(int i=1; i<argc; i++)
		args.push_back(argv[i]);
#ifdef WITH_OPENGL
	vsr::initWindow("RotationShield Visualizer");
	pthread_t thread;
	pthread_create( &thread, NULL, &menuThread, &args );
	vsr::doGlutLoop();
#else
	menuSystem(args);
#endif
	return 0;
}
