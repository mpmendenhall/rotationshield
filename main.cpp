/// file "main.cpp" \brief Example of shield simulation


#include "Sandbox.hh"
#include "tests.hh"
#include <iostream>
#include <fstream>
#include "ncube.hh"
#include "VarMat.hh"

#include "HalfScaleCoil.hh"
#include "FullScaleCoil.hh"

int main (int argc, char * const argv[]) {
	
	cubeDemo();
	exit(0);
	
	std::string basepath = "/Users/michael/Documents/EDM/ShieldStudies/";
	
	std::string projectpath = basepath+"/Bare_Fullscale/";
	makeDir(projectpath);
	
	std::string fieldspath = basepath+"/Fields/";
	makeDir(fieldspath);
	
	std::string gfpath = basepath+"/GFs/";
	makeDir(gfpath);
	
	std::ofstream fieldsout;
	std::ofstream statsout;
	fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
	statsout.open((fieldspath+"/Fieldstats.txt").c_str());
	
	bareFullScaleCoil(fieldsout, statsout);
	//fullScaleCoil(fieldsout, statsout);
	//halfScaleCoil(0.03175,fieldsout,statsout);
	
	fieldsout.close();
	statsout.close();
	
	vsr::Visr::W->pause();
	
	return 0;
}





int main_halfstudy (int argc, char * const argv[]) {
	
	const char* basepath = "/Users/michael/Documents/EDM/ShieldStudies/";
	
	char projectpath[1024];
	sprintf(projectpath,"%s/Halfscale_TiltY/",basepath);
	makeDir(projectpath);
	
	char fieldspath[1024];
	sprintf(fieldspath,"%s/Fields/",projectpath);
	makeDir(fieldspath);
	
	char gfpath[1024];
	sprintf(gfpath,"%s/GFs/",projectpath);
	makeDir(gfpath);
	
	std::ofstream fieldsout;
	std::ofstream statsout;
	
	char tmpc[1024];
	
	ScanRange sr(-0.01,0.01,11);
	
	for(mdouble a = sr.next(); sr.goOn(); a=sr.next()) {
		
		printf("Running calculation for dr = "); sr.printStatus(); printf("\n");
		
		sprintf(tmpc,"%s/Fieldmap_%f.txt",fieldspath,(double)a);
		fieldsout.open(tmpc);
		sprintf(tmpc,"%s/Fieldstats_%f.txt",fieldspath,(double)a);
		statsout.open(tmpc);

		//halfScaleTranslated(vec3(0,a,0), fieldsout, statsout);
		halfScaleTilted(vec3(0,a,0), fieldsout, statsout);
		//vsr::Visr::W->pause();
		
		fieldsout.close();
		statsout.close();
	}
	
	vsr::Visr::W->pause();
	
	return 0;
}

