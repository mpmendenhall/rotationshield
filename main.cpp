/// file "main.cpp" \brief Menu system for launching various shield simulations

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

/// Return a random number, uniformly distributed over interval [a,b]
/**	\param a lower bound of interval
 \param b upper bound of interval
 \return a random number in the interval [a,b] */
double randunif(double a, double b) { 
	return a + (b-a)*double(rand())/double(RAND_MAX);
}

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
	printf("Simple shield self-test...\n");
	reference_simpleshield();
}

/// run demonstration routines
void mi_demos(StreamInteractor* S) {
	std::string s = S->popString();
	if(s=="sc")
		superball_test();
	if(s=="mr")
		mirror_test();
	if(s=="tb")
		tube_test();
	if(s=="ft")
		flux_trap_test();
}

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
	
	InputRequester exitMenu("Exit Menu", &menutils_Exit);
	InputRequester ncube("Hyercube visualization test", &mi_ncube);
	ncube.addArg("n dim","3");
	InputRequester setClearColor("Set visualization background color", &mi_clearcolor);
	setClearColor.addArg("r","0.0");
	setClearColor.addArg("g","0.0");
	setClearColor.addArg("b","0.0");
	
	InputRequester selfTests("Self test to reproduce known results", &mi_runtests);
	
	NameSelector selectDemo("Demonstration");
	selectDemo.addChoice("superconducting sphere","sc");
	selectDemo.addChoice("superconductor-mirrored cos theta coil","mr");
	selectDemo.addChoice("cos-theta coil in ferromagnetic tube","tb");
#ifdef WITH_LAPACKE
	selectDemo.addChoice("trapped-flux state of superconducting ring","ft");
#endif
	InputRequester demos("Demonstration calculations", &mi_demos);
	demos.addArg(&selectDemo);
	
	
	OptionsMenu OM("Rotation Shield Main Menu");
#ifdef WITH_OPENGL
	OM.addChoice(&ncube,"ncube",SELECTOR_HIDDEN);
	OM.addChoice(&setClearColor,"cc",SELECTOR_HIDDEN);
#endif
	OM.addChoice(&demos,"demo");
	
	SystemConfiguration SC;
	OM.addChoice(&SC.outDir,"dir");
	OM.addChoice(&SC.OMfieldsrc,"field");
	OM.addChoice(&SC.OMsurfaces,"bound");
	OM.addChoice(&SC.OMcell,"cell");
	
	OM.addChoice(&SC.doSolve,"solve");
	OM.addChoice(&SC.doApply,"apply");
	OM.addChoice(&SC.doMeas,"meas");
#ifdef WITH_LAPACKE
	OM.addChoice(&SC.addSingular,"svd");
	OM.addChoice(&SC.setSingularEpsilon,"ep");
#endif
	
	OM.addChoice(&exitMenu,"x");
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
