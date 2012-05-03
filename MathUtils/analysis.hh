/// \file "analysis.hh" \brief Routines for analyzing the simulated fields
#ifndef ANALYSIS_HH
/// Make sure this header is only loaded once
#define ANALYSIS_HH 1

#include "FieldSource.hh"
#include "VarVec.hh"
#include "linmin.hh"
#include "Geometry.hh"
#include <iostream>
#include <iomanip>

/// A class for analyzing field uniformity/shape over rectangular grids.
///
/// Uses Simpson's Rule to numerically integrate the fields
///	over the specified volume --- good enough for the smoothly
///	varying fields typically encountered.
class FieldAnalyzer {
public:
	/// Constructor \param fs the FieldSource to analyze fields from
	FieldAnalyzer(FieldSource* fs): FS(fs) { FS->retain(); }
	/// destructor
	~FieldAnalyzer() { FS->release(); }
	/// Survey the fields over a rectangular grid
	/// \param ll lower-left corner of the survey area
	/// \param ur epur-right corner of the survey area
	/// \param nX number of grid points along the x direction
	/// \param nY number of grid points along the y direction
	/// \param nZ number of grid points along the z direction
	/// \param datf file to store surveyed fields to
	/// \param statsout ostream to which summary fields data is written
	/// \param key string to prepend to summary data
	void survey(vec3 ll, vec3 ur, int nX, int nY, int nZ, std::ostream& statsout = std::cout, std::ostream& datsout = std::cout) const;
	
	void visualizeSurvey(vec3 ll, vec3 ur, int nX, int nY, int nZ) const;
	
	/// Get the Simpson's Rule numerical integrating weight for a point
	static mdouble simpsonCoeff(int n, int ntot);
	/// Get the Simpson's Rule numerical integrating weight for a point
	static mdouble simpsonCoeff(mdouble n, int ntot) { return simpsonCoeff(int(n),ntot); }	
protected:
	const FieldSource* FS; //< the FieldSource producing the fields being analyzed
};

mdouble FieldAnalyzer::simpsonCoeff(int n, int ntot) {
	if(ntot == 1) return 1.0;
	if(ntot == 2) return 0.5;
	if(n==0 || n==ntot-1) return 1.0/(3.0*mdouble(ntot-1));
	if(n%2==0) return 2.0/(3.0*mdouble(ntot-1));
	return 4.0/(3.0*mdouble(ntot-1));
}

void FieldAnalyzer::survey(vec3 ll, vec3 ur, int nX, int nY, int nZ, std::ostream& statsout, std::ostream& datsout) const {	
	printf("Analyzing resultant fields... ");
	vec3 n = vec3(nX-1,nY-1,nZ-1);
	vec3 dx = (ur - ll)/n;
	for(int i=0; i<3; i++) if(n[i]==0) { dx[i] = 0; ll[i] = 0.5*(ll[i]+ur[i]); }
	vec3 b, x, bsZ, bsYZ, bsXYZ, bssZ, bssYZ, bssXYZ;
	mdouble avgBx0 = 0;
	mdouble avgBx1 = 0;
	
	// set up results tables
	CLHEP::HepMatrix coords(nX*nY*nZ,3);
	CLHEP::HepVector bfield[3];
	for(int i=0; i<3; i++) bfield[i] = CLHEP::HepVector(nX*nY*nZ);
	
	// scan over survey points
	bsXYZ = vec3(); bssXYZ = vec3();
	int c=0;
	for(n[0] = 0; n[0] < nX; n[0]++) {
		bsYZ = vec3(); bssYZ = vec3();
		for(n[1] = 0; n[1] < nY; n[1] ++) {
			
			printf("*"); fflush(stdout);
			
			bsZ = vec3(); bssZ = vec3();
			for(n[2] = 0; n[2] < nZ; n[2]++)
			{
				x = ll+n*dx;
				b = FS->fieldAt(x);
				for(int i=0; i<3; i++) bfield[i][c] = b[i];
				
				bsZ += b*simpsonCoeff(n[2],nZ);
				bssZ += b*b*simpsonCoeff(n[2],nZ);
				
				for(int i=0; i<3; i++) datsout << std::setw(16) << x[i] << " ";
				for(int i=0; i<3; i++) datsout << std::setw(20) << std::setprecision(12) << bfield[i][c] << " ";
				datsout << std::endl;
				
				for(int i=0; i<3; i++) coords[c][i] = x[i];
				
				c++;
			}
			bsYZ += bsZ*simpsonCoeff(n[1],nY);
			bssYZ += bssZ*simpsonCoeff(n[1],nY);
		}
		if(n[0]==0) avgBx0 = bsYZ[0];
		if(n[0]==nX-1) avgBx1 = bsYZ[0];
		bsXYZ += bsYZ*simpsonCoeff(n[0],nX);
		bssXYZ += bssYZ*simpsonCoeff(n[0],nX);
	}
	printf("\n");
	
	
	// summary data and field fit
	statsout << "#Field Shape" << std::endl;
	statsout << "#<B> = (" << bsXYZ[0] << " " << bsXYZ[1] << " " << bsXYZ[2] << ")\n";
	statsout << "#<B^2> = (" << bssXYZ[0] << " " << bssXYZ[1] << " " << bssXYZ[2] << ")\n";
	statsout << "#sigmaBx/Bx0 = " << sqrt(bssXYZ[0] - bsXYZ[0]*bsXYZ[0])/bsXYZ[0] << std::endl;
	statsout << "#<dBx/dx>/Bx0 = " << (avgBx1-avgBx0)/(bsXYZ[0]*(ur[0]-ll[0])) << std::endl;
	mdouble resids;
	char axisnames[] = "xyz";
	
	Polynomial<3,mdouble> p = Polynomial<3,mdouble>::lowerTriangleTerms(4);
	for(int i=0; i<3; i++) {
		resids = polynomialFit(coords, bfield[i], p);
		statsout << "#4th order polynomial fit to B" << axisnames[i] << " centered at < 0 0 0 >, rms residual = " << resids << std::endl;
		Polynomial<3,mdouble> p2 = p;
		p2.prune(1e-9*fabs(p(vec3())));
		p2.tableForm(statsout);
		//statsout << "#4th order Polynomial for B" << axisnames[i] << " recentered at " << (ll+ur)*0.5 << std::endl;
		//p2 = p.recentered((ll+ur)*0.5);
		//p2.prune(1e-9);
		//p2.tableForm(statsout);
	}
	
}


void FieldAnalyzer::visualizeSurvey(vec3 ll, vec3 ur, int nX, int nY, int nZ) const {	
	
	vec3 n = vec3(nX-1,nY-1,nZ-1);
	vec3 dx = (ur - ll)/n;
	vec3 x,b;
	for(int i=0; i<3; i++) if(n[i]==0) { dx[i] = 0; ll[i] = 0.5*(ll[i]+ur[i]); }
	
	vsr::Visr::W->startRecording(); 
	vsr::Visr::W->clearWindow();
	FS->visualize(false);
	
	for(n[0] = 0; n[0] < nX; n[0]++) {
		for(n[1] = 0; n[1] < nY; n[1] ++) {
			for(n[2] = 0; n[2] < nZ; n[2]++) {
				
				x = ll+n*dx;
				b = FS->fieldAt(x);
				float bmag = b.mag();
				if(!bmag) continue;
				b *= log(1+pow(bmag,0.25))/bmag;
				Line(x-b*0.02,x+b*0.02).visualizeDirected();
			}
		}
	}
	
	vsr::Visr::W->stopRecording();
}


#endif
