/// file "main.cpp" \brief Menu system for launching various shield simulations

#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>

#include "PathUtils.hh"
#include "MiscUtils.hh"
#include "ControlMenu.hh"
#include "tests.hh"
#include "Studies.hh"
#include "ncube.hh"

/// run self-tests
void mi_runtests(std::deque<std::string>&, std::stack<std::string>&) {
	printf("Simple shield self-test...\n");
	reference_simpleshield();
	printf("\n\nReference self-test...\n");
	reference_sanity_check();
}

/// hypercube visualization test
void mi_ncube(std::deque<std::string>&, std::stack<std::string>&) {
	const unsigned int N = 4;
	NRotator<N> NR;
	NCube<N> NC;
	
	double ftime = 15.e-3; // frame time in s
	std::vector<mdouble> vrot;
	for(unsigned int i=1; i<N; i++)
		for(unsigned int j=0; j<i; j++)
			vrot.push_back(randunif(-1.*ftime,1.*ftime));
	
	vsr::set_pause();
	if(vsr::get_pause())
		printf("Press [ENTER] in visualization window to continue...\n");
	while(vsr::get_pause()) {
		unsigned int n=0;
		for(unsigned int i=1; i<N; i++)
			for(unsigned int j=0; j<i; j++)
				NR.rotate(i, j, vrot[n++]);
		NC.visualize(&NR);
		usleep(ftime*1.e6);
	}
}

// Global coil/shield settings
nEDM_Geom* GlobGM = NULL;
CosThetaBuilder* GlobCT = NULL;
coilShield* GlobCS = NULL;
fieldCell* GlobFC = NULL;

// set cell measurement range
void mi_setFCrange(std::deque<std::string>&, std::stack<std::string>& stack) {
	float urz = streamInteractor::popFloat(stack);
	float ury = streamInteractor::popFloat(stack);
	float urx = streamInteractor::popFloat(stack);
	float llz = streamInteractor::popFloat(stack);
	float lly = streamInteractor::popFloat(stack);
	float llx = streamInteractor::popFloat(stack);
	if(!GlobFC) GlobFC = new fieldCell(vec3(),vec3(),5,5,5);
	GlobFC->ll = vec3(llx,lly,llz);
	GlobFC->ur = vec3(urx,ury,urz);
	GlobFC->getInfo().display();
}
// set cell gridding
void mi_setFCgrid(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobFC) GlobFC = new fieldCell(vec3(),vec3(),5,5,5);
	GlobFC->nz = streamInteractor::popInt(stack);
	GlobFC->ny = streamInteractor::popInt(stack);
	GlobFC->nx = streamInteractor::popInt(stack);
	GlobFC->getInfo().display();
}

// set coil parameters
void mi_setCoil(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobCT) GlobCT = new CosThetaBuilder(0,0,0);
	GlobCT->radius = streamInteractor::popFloat(stack);
	GlobCT->length = streamInteractor::popFloat(stack);
	GlobCT->ncoils = streamInteractor::popInt(stack);
}
// set coil distortion
void mi_setCoilDistort(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobCT) GlobCT = new CosThetaBuilder(0,0,0);
	float a = streamInteractor::popFloat(stack);
	int n = streamInteractor::popInt(stack);
	ShiftPositioner* SP = (ShiftPositioner*)(GlobCT->AP);
	if(n<=0) SP->shift = VarVec<mdouble>(0);
	else {
		while(SP->shift.size()<(unsigned int)n) SP->shift.push_back(0);
		SP->shift[n-1] = a;
	}
}

// set shield parameters
void mi_setShield(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobGM) GlobGM = new nEDM_Geom();
	if(!GlobCS) GlobCS = new coilShield();
	GlobGM->shield = GlobCS;
	GlobCS->radius = streamInteractor::popFloat(stack);
	GlobCS->length = streamInteractor::popFloat(stack);
}
// set unshielded
void mi_unshield(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(GlobGM) GlobGM->shield = NULL;
}
// set shield gridding
void mi_shieldGrid(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobCS) GlobCS = new coilShield();
	GlobCS->pSegs = streamInteractor::popInt(stack);
	GlobCS->vSegs = streamInteractor::popInt(stack);
	GlobCS->cSegs = streamInteractor::popInt(stack);
}

// set save grid
void mi_setSaveGrid(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobGM) GlobGM = new nEDM_Geom();
	GlobGM->saveGrid = streamInteractor::popInt(stack);
}

// take measurement
void mi_meas(std::deque<std::string>&, std::stack<std::string>& stack) {
	if(!GlobGM) { GlobGM = new nEDM_Geom(); }
	GlobGM->basedir = getEnvSafe("ROTSHIELD_OUT","")+"/"+streamInteractor::popString(stack);
	if(!GlobCT) {
		printf("Coil not specified! Nothing to measure!\n");
		return;
	}
	if(!GlobFC) {
		printf("Cell not specified! Nothing to measure!\n");
		return;
	}
	if(!GlobGM->shield) printf("Measuring bare coil\n");
	else printf("Measuring shielded coil\n");
	GlobGM->coil = GlobCT;
	GlobGM->construct();
	GlobGM->cell = GlobFC;
	GlobGM->takeSample();
	printf("Data collection complete.\n");
}

/// Short ("vertical") coil with re-oriented sample cell
/// N=30, length=2.m, radius variable
/// shield radius nominal 10cm outside coil; shield 40cm longer than coil
/// cell dimensions y=40cm, x=7.5cm, z=10cm, inner edge 5cm from center in x (field) direction
/// Want:
///   B0 = 30mG; <dBi/di> < 10^-7 G/cm
///   sqrt(<(dBx/dz)^2>) (z=long direction, x=along field) small
///   cells: | 7.5 | 10 | 7.5 |    --> x


void menuSystem(std::deque<std::string> args=std::deque<std::string>()) {

	GlobFC = new fieldCell(vec3(-.2,-.2,-.2),vec3(.2,.2,.2),11,11,11);
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester peek("Show stack",&menutils_PrintStack);
	inputRequester ncube("Hyercube visualization test",&mi_ncube);
	inputRequester selfTests("Self test to reproduce known results",&mi_runtests);
	
	inputRequester setFCrange("Set Measurement Range",&mi_setFCrange);
	setFCrange.addArg("x min","-0.20");
	setFCrange.addArg("y min","-0.20");
	setFCrange.addArg("z min","-0.20");
	setFCrange.addArg("x max","0.20");
	setFCrange.addArg("y max","0.20");
	setFCrange.addArg("z max","0.20");
	inputRequester setFCgrid("Set Measurement Gridding",&mi_setFCgrid);
	setFCgrid.addArg("nx","5");
	setFCgrid.addArg("ny","5");
	setFCgrid.addArg("nz","5");
	OptionsMenu OMcell("Cell Options");
	OMcell.addChoice(&setFCrange,"range");
	OMcell.addChoice(&setFCgrid,"grid");
	OMcell.addChoice(&exitMenu,"x");
	
	inputRequester setCoil("Set Coil Geometry",&mi_setCoil);
	setCoil.addArg("half n","15");
	setCoil.addArg("Length","4.292");
	setCoil.addArg("Radius","0.61");
	inputRequester setCoilDistort("Set Coil Distortion Parameter",&mi_setCoilDistort);
	setCoilDistort.addArg("param","0");
	setCoilDistort.addArg("value","0");
	OptionsMenu OMcoil("Coil Options");
	OMcoil.addChoice(&setCoil,"geom");
	OMcoil.addChoice(&setCoilDistort,"dist");
	OMcoil.addChoice(&exitMenu,"x");
	
	inputRequester setShield("Set Shield Geometry",&mi_setShield);
	setShield.addArg("Length","4.692");
	setShield.addArg("Radius","0.68");
	inputRequester unsetShield("Remove shield",&mi_unshield);
	inputRequester setShieldGrid("Set Shield Gridding",&mi_shieldGrid);
	setShieldGrid.addArg("Fixed Z segs","10");
	setShieldGrid.addArg("Variable Z segs","20");
	setShieldGrid.addArg("Phi segs","128");
	OptionsMenu OMshield("Shield Options");
	OMshield.addChoice(&setShield,"geom");
	OMshield.addChoice(&setShieldGrid,"grid");
	OMshield.addChoice(&unsetShield,"clear");
	OMshield.addChoice(&exitMenu,"x");
	
	inputRequester meas("Coil field measurement",&mi_meas);
	meas.addArg("Output Directory");
	inputRequester setSaveGrid("Enable/disable gridded output",&mi_setSaveGrid);
	setSaveGrid.addArg("Enable","1");
	OptionsMenu OMmeas("Measurement Routines");
	OMmeas.addChoice(&meas,"run");
	OMmeas.addChoice(&setSaveGrid,"svgrd");
	OMmeas.addChoice(&exitMenu,"x");
	
	OptionsMenu OM("Rotation Shield Main Menu");
	//OM.addChoice(&selfTests,"test");
#ifdef WITH_OPENGL
	OM.addChoice(&ncube,"ncube");
#endif
	OM.addChoice(&OMcell,"cell");
	OM.addChoice(&OMcoil,"coil");
	OM.addChoice(&OMshield,"shield");
	OM.addChoice(&OMmeas,"meas");
	OM.addChoice(&exitMenu,"x");
	OM.addChoice(&peek,"peek",SELECTOR_HIDDEN);
	std::stack<std::string> stack;
	OM.doIt(args,stack);
	
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
	vsr::initWindow();
	pthread_t thread;
	pthread_create( &thread, NULL, &menuThread, &args );
	vsr::doGlutLoop();
#else
	menuSystem(args);
#endif
	return 0;
}
