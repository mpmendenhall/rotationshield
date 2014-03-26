/// \file Studies.hh \brief Utility classes for setting up Cos Theta Coil studies

#ifndef STUDIES_HH
/// Makes sure to only load this file once
#define STUDIES_HH 1

#include "MixedSource.hh"
#include "SymmetricSolver.hh"
#include "CosThetaBuilder.hh"
#include "MagRS.hh"
#include "ControlMenu.hh"
#include "QFile.hh"
#include "DVFunc.hh"
#include <string>

/// Data-sampling cell specification
class fieldCell {
public:
	/// default constructor
	fieldCell(): ll(-.2,-.2,-.2), ur(.2,.2,.2), nx(5), ny(5), nz(5), vx(2), vy(3), vz(2), saveGrid(true) {}
	/// cell info as Stringmap
	Stringmap getInfo() const;
	
	vec3 ll;				//< lower corner
	vec3 ur;				//< upper corner
	unsigned int nx,ny,nz;	//< number of points to sample along each axis
	unsigned int vx,vy,vz;	//< number of points to visualize along each axis
	bool saveGrid;			//< whether to save gridded output on survey
};

/// Utility class for configuring field measurement
class SystemConfiguration: public StreamInteractor {
public:
	/// constructor
	SystemConfiguration();
	/// destructor
	~SystemConfiguration();
	
	/// initialize symmetric boundary system
	void initReactiveSet(unsigned int nPhi);
	
	/// adapt a geometry to current field source
	DVFunc1<2,double>* adaptSurface(DVFunc1<2,double>* f, double pfixed, bool useTotal = false) const;
	/// solve surface interaction system
	void solve(const std::string& cfile = "");
	/// calculate applied response
	void calculate_result();
	
	/// measure fields, writing to given output file
	void measureFields(const std::string& xpath="") const;
	/// write measurement info to file at basedir/xpath/GeomInfo.txt
	void writeInfo(const std::string& xpath="") const;
	
	std::string basedir;			//< base directory for IO operations
	MagRSCombiner* RSC;				//< reacting boundary condition surfaces
	MixedSource* IncidentSource;	//< incident field source
	MixedSource* TotalField;		//< incident + reacting field
	fieldCell cell;					//< field measurement cell
	SymmetricSolver* SS;			//< system solver
	CosThetaBuilder CTB;			//< cos theta field coil builder
	
	
	
	//------------------------------------
	// menu-driven user interface elements
	
	InputRequester exitMenu;
	InputRequester outDir;
	
	InputRequester setFCrange;
	InputRequester setFCgrid;
	InputRequester setSaveGrid;
	OptionsMenu OMcell;
	
	InputRequester addLineCurrent;
	InputRequester addUnifB;
	InputRequester clearIncident;
	InputRequester buildCosThetaExit;
	OptionsMenu OMfieldsrc;
	
	InputRequester setPhi;
	InputRequester addSlab;
	InputRequester addTube;
	//InputRequester addSheet;
	OptionsMenu OMsurfaces;
	
	InputRequester doSolve;
	InputRequester doApply;
	InputRequester doMeas;
};

#endif
