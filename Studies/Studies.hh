/// \file Studies.hh \brief Utility classes for setting up Cos Theta Coil studies

#ifndef STUDIES_HH
/// Makes sure to only load this file once
#define STUDIES_HH 1

#include "CosThetaBuilder.hh"
#include "SymmetricSolver.hh"
#include "PlaneSource.hh"
#include "ShieldBuilder.hh"
#include "analysis.hh"
#include "QFile.hh"

/// Generic shield for cos theta coil
class coilShield {
public:
	/// constructor
	coilShield(): length(0), radius(0), mu(10000.), cSegs(10), vSegs(20), pSegs(128) {}
	/// shield info as Stringmap
	Stringmap getInfo() const;
	/// construct into mixed source
	void construct(MixedSource& ms, CosThetaBuilder* ct) const;
	
	mdouble length;			//< shield length
	mdouble radius;			//< shield radius
	double mu;				//< material relative permeability
	unsigned int cSegs;		//< constant-width z segments
	unsigned int vSegs;		//< variable z segments
	unsigned int pSegs;		//< phi segments
};

/// Data-sampling cell
class fieldCell {
public:
	/// constructor
	fieldCell(vec3 l, vec3 u, unsigned int n1, unsigned int n2, unsigned int n3): ll(l), ur(u), nx(n1), ny(n2), nz(n3), vx(3), vy(3), vz(3) {}
	/// cell info as Stringmap
	Stringmap getInfo() const;
	
	vec3 ll;				//< lower corner
	vec3 ur;				//< upper corner
	unsigned int nx,ny,nz;	//< number of points to sample along each axis
	unsigned int vx,vy,vz;	//< number of points to visualize along each axis
};

/// Generic nEDM cos-theta coil and shield geometry measurement
class nEDM_Geom {
public:
	/// constructor
	nEDM_Geom(const std::string& dir = "");
	/// destructor
	virtual ~nEDM_Geom();
	/// construct coil/shield and calculate
	void construct();
	/// sample fields
	void takeSample(const std::string& sName="Fields");
	
	CosThetaBuilder* coil;
	coilShield* shield;
	fieldCell* cell;
	
	bool saveGrid;			//< whether to save the grid data to file
	std::string basedir;	//< directory for output
	
protected:
	MixedSource* ms;		//< combined coil+shield field source
	
};

#endif
