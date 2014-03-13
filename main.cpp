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

/// run self-tests
void mi_runtests(std::deque<std::string>&, std::stack<std::string>&) {
	csurface_test();
	return;
	printf("Simple shield self-test...\n");
	reference_simpleshield();
	//printf("\n\nReference self-test...\n");
	//reference_sanity_check();
}

/// Return a random number, uniformly distributed over interval [a,b]
/**	\param a lower bound of interval
 \param b upper bound of interval
 \return a random number in the interval [a,b] */
mdouble randunif(mdouble a, mdouble b) { 
	return a + (b-a)*mdouble(rand())/mdouble(RAND_MAX);
}

template<unsigned int N>
void hypercube_rotator() {
	NRotator<N> NR;
	NCube<N> NC;
	
	double ftime = 15.e-3; // min frame time in s
	srand(time(NULL));
	std::vector<mdouble> vrot;
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

/// hypercube visualization test
void mi_ncube(std::deque<std::string>&, std::stack<std::string>& stack) {
	int n = streamInteractor::popInt(stack);
	if(n==3) hypercube_rotator<3>();
	else if(n==4) hypercube_rotator<4>();
	else if(n==5) hypercube_rotator<5>();
	else printf("Only 3 to 5 dimensions allowed.\n");
}

void mi_clearcolor(std::deque<std::string>&, std::stack<std::string>& stack) {
	float b = streamInteractor::popFloat(stack);
	float g = streamInteractor::popFloat(stack);
	float r = streamInteractor::popFloat(stack);
	vsr::setClearColor(r,g,b);
	vsr::startRecording();
	vsr::clearWindow();
	vsr::stopRecording();
}



// Global coil/shield settings
nEDM_Geom GlobGM;

// set cell measurement range
void mi_setFCrange(std::deque<std::string>&, std::stack<std::string>& stack) {
	float urz = streamInteractor::popFloat(stack);
	float ury = streamInteractor::popFloat(stack);
	float urx = streamInteractor::popFloat(stack);
	float llz = streamInteractor::popFloat(stack);
	float lly = streamInteractor::popFloat(stack);
	float llx = streamInteractor::popFloat(stack);
	GlobGM.cell.ll = vec3(llx,lly,llz);
	GlobGM.cell.ur = vec3(urx,ury,urz);
}
// set cell gridding
void mi_setFCgrid(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.cell.nz = streamInteractor::popInt(stack);
	GlobGM.cell.ny = streamInteractor::popInt(stack);
	GlobGM.cell.nx = streamInteractor::popInt(stack);
}

// set coil parameters
void mi_setCoil(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.coil.radius = streamInteractor::popFloat(stack);
	GlobGM.coil.length = streamInteractor::popFloat(stack);
	GlobGM.coil.ncoils = streamInteractor::popInt(stack);
}
// set coil distortion
void mi_setCoilDistort(std::deque<std::string>&, std::stack<std::string>& stack) {
	float a = streamInteractor::popFloat(stack);
	int n = streamInteractor::popInt(stack);
	ShiftPositioner* SP = (ShiftPositioner*)(GlobGM.coil.AP);
	if(n<=0) SP->shift = VarVec<mdouble>(0);
	else {
		while(SP->shift.size()<(unsigned int)n) SP->shift.push_back(0);
		SP->shift[n-1] = a;
	}
}
// set coil end
void mi_setCoilEnd(std::deque<std::string>&, std::stack<std::string>& stack) {
	std::string side = streamInteractor::popString(stack);
	std::string endtp = streamInteractor::popString(stack);
	for(unsigned int zside = 0; zside<2; zside++) {
		if((side=="neg" && zside==0) || (side=="pos" && zside==1)) continue;
		if(endtp=="arc") GlobGM.coil.myCap[zside] = CosThetaBuilder::CAP_ARC;
		if(endtp=="line") GlobGM.coil.myCap[zside] = CosThetaBuilder::CAP_LINE;
		if(endtp=="none") GlobGM.coil.myCap[zside] = CosThetaBuilder::CAP_NONE;
	}
}

// set shield frame parameters
void mi_setShield(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.shield.pSegs = streamInteractor::popInt(stack);
	GlobGM.shield.radius = streamInteractor::popFloat(stack);
	GlobGM.shield.length = streamInteractor::popFloat(stack);
}

GeomRefPt refPtByName(const std::string& s) {
	if(s=="nax") return GEOMREF_NEGAXIS;
	if(s=="pax") return GEOMREF_POSAXIS;
	if(s=="nrd") return GEOMREF_NEGRADIUS;
	if(s=="prd") return GEOMREF_POSRADIUS;
	if(s=="c") return GEOMREF_CENTER;
	return GEOMREF_ORIGIN;
}

// build new shield section
void mi_newShieldSection(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.shield.mySections.push_back(shieldSection());
	GlobGM.shield.mySections.back().endpts[1] = refPtByName(streamInteractor::popString(stack));
	GlobGM.shield.mySections.back().endpts[0] = refPtByName(streamInteractor::popString(stack));
	GlobGM.shield.mySections.back().vSegs = streamInteractor::popInt(stack);
	GlobGM.shield.mySections.back().cSegs = streamInteractor::popInt(stack);
	GlobGM.shield.mySections.back().mu = streamInteractor::popFloat(stack);
}

// modify section offset
void mi_shieldOffsets(std::deque<std::string>&, std::stack<std::string>& stack) {
	assert(GlobGM.shield.mySections.size());
	GlobGM.shield.mySections.back().endoff[1][1] = streamInteractor::popFloat(stack);
	GlobGM.shield.mySections.back().endoff[1][0] = streamInteractor::popFloat(stack);
	GlobGM.shield.mySections.back().endoff[0][1] = streamInteractor::popFloat(stack);
	GlobGM.shield.mySections.back().endoff[0][0] = streamInteractor::popFloat(stack);
}

// set save grid
void mi_setSaveGrid(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.saveGrid = streamInteractor::popInt(stack);
}

// take measurement
void mi_meas(std::deque<std::string>&, std::stack<std::string>& stack) {
	GlobGM.basedir = getEnvSafe("ROTSHIELD_OUT","")+"/"+streamInteractor::popString(stack);
	if(!GlobGM.coil.ncoils) {
		printf("Coil not specified! Nothing to measure!\n");
		return;
	}
	if(!GlobGM.shield.mySections.size()) printf("Measuring bare coil\n");
	else printf("Measuring shielded coil\n");
	GlobGM.construct();
	GlobGM.takeSample();
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
	
	inputRequester exitMenu("Exit Menu",&menutils_Exit);
	inputRequester ncube("Hyercube visualization test",&mi_ncube);
	ncube.addArg("n dim","3");
	inputRequester setClearColor("Set visualization background color",&mi_clearcolor);
	setClearColor.addArg("r","0.0");
	setClearColor.addArg("g","0.0");
	setClearColor.addArg("b","0.0");
	
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
	NameSelector selectCoilEnd("Coil end wire shape");
	selectCoilEnd.addChoice("smooth arc","arc");
	selectCoilEnd.addChoice("straight line","line");
	selectCoilEnd.addChoice("no end wires","none");
	selectCoilEnd.setDefault("arc");
	NameSelector selectSide("z side");
	selectSide.addChoice("negative","neg");
	selectSide.addChoice("positive","pos");
	selectSide.addChoice("both","both");
	selectSide.setDefault("both");
	inputRequester setCoilEnd("End wires",&mi_setCoilEnd);
	setCoilEnd.addArg(&selectCoilEnd);
	setCoilEnd.addArg(&selectSide);
	OptionsMenu OMcoil("Coil Options");
	OMcoil.addChoice(&setCoil,"geom");
	OMcoil.addChoice(&setCoilDistort,"dist");
	OMcoil.addChoice(&setCoilEnd,"end");
	OMcoil.addChoice(&exitMenu,"x");
	
	inputRequester setShield("Set frame geometry",&mi_setShield);
	setShield.addArg("length");
	setShield.addArg("radius");
	setShield.addArg("phi grid segments","128");
	
	NameSelector selectRefPt("Start reference point");
	selectRefPt.addChoice("negative axis","nax");
	selectRefPt.addChoice("positive axis","pax");
	selectRefPt.addChoice("negative radius","nrd");
	selectRefPt.addChoice("positive radius","prd");
	selectRefPt.addChoice("center","c");
	selectRefPt.setDefault("nrd");
	NameSelector selectEndRefPt = selectRefPt;
	selectEndRefPt.name = "End reference point";
	selectEndRefPt.setDefault("prd");
	inputRequester addShieldSect("Add section",&mi_newShieldSection);
	addShieldSect.addArg("mu","10000");
	addShieldSect.addArg("Constant Z segs","10");
	addShieldSect.addArg("Variable Z segs","20");
	addShieldSect.addArg(&selectRefPt);
	addShieldSect.addArg(&selectEndRefPt);
	
	inputRequester modShieldOffset("Modify section geometry",&mi_shieldOffsets);
	modShieldOffset.addArg("start dz","0");
	modShieldOffset.addArg("start dr","0");
	modShieldOffset.addArg("end dz","0");
	modShieldOffset.addArg("end dr","0");
	
	OptionsMenu OMshield("Shield Options");
	OMshield.addChoice(&setShield,"geom");
	OMshield.addChoice(&addShieldSect,"add");
	OMshield.addChoice(&modShieldOffset,"mod");
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
#ifdef WITH_OPENGL
	OM.addChoice(&ncube,"ncube");
	OM.addChoice(&setClearColor,"bg");
#endif
	OM.addChoice(&OMcell,"cell");
	OM.addChoice(&OMcoil,"coil");
	OM.addChoice(&OMshield,"shield");
	OM.addChoice(&OMmeas,"meas");
	OM.addChoice(&exitMenu,"x");
	OM.addChoice(&selfTests,"test",SELECTOR_HIDDEN);
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
	vsr::initWindow("RotationShield Visualizer");
	pthread_t thread;
	pthread_create( &thread, NULL, &menuThread, &args );
	vsr::doGlutLoop();
#else
	menuSystem(args);
#endif
	return 0;
}
