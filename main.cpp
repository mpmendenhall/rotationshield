/// file "main.cpp" \brief Example of shield simulation

#include <iostream>
#include <fstream>

#include "tests.hh"
#include "HalfScaleCoil.hh"
#include "FullScaleCoil.hh"

int main (int argc, char * const argv[]) {
	
	reference_sanity_check();
	
	// set up output paths
	std::string basepath = "/Users/michael/Documents/EDM/ShieldStudies/";
	std::string projectpath = basepath+"/Bare_Fullscale/";
	makeDir(projectpath);
	std::string fieldspath = basepath+"/Fields/";
	makeDir(fieldspath);
	std::string gfpath = basepath+"/GFs/";
	makeDir(gfpath);
	
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
	
	// optional: construct shield
	if(0) {
		mdouble shieldDistance = 0.069;
		mdouble sradius = cradius + shieldDistance;
		mdouble slen = clen + 0.40;
		FieldEstimator2D* fe = new FieldEstimator2D();
		fe->addsource(vec2(-clen/2.0,cradius),1.0);
		fe->addsource(vec2(clen/2.0,cradius),1.0);
		ShieldBuilder* G = new ShieldBuilder(128);
		G->makeOptCyl(10, 20, sradius, -slen/2, slen/2, new PlaneSource(Plane(),10000), fe);
		SymmetricSolver* sp = new SymmetricSolver(G);
		sp->solve();
		sp->calculateIncident(ms);
		sp->calculateResult();
		ms->addsource(sp);
		ms->visualize();
	}

	
	// survey data points on grid
	FieldAnalyzer FA = FieldAnalyzer(ms);
	//        lower left corner,     upper right corner, nx,ny,nz, output files
	FA.survey(vec3(0.0,0.0,0.0),vec3(0.20,0.20,0.5),10,10,10,statsout,fieldsout);
	
	
	// cleanup
	ms->release();
	fieldsout.close();
	statsout.close();
	vsr::Visr::W->pause();
	
	return 0;
}
