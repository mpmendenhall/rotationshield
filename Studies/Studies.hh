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
	coilShield();
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
	
	// optional endcap on each side
	mdouble endcap_ir[2];			//< inner radius for (optional) superconducting endcap
	mdouble endcap_delta_or[2];		//< difference between endcap outer radius and main shield
	mdouble endcap_delta_z[2];		//< endcap offset from end of main shield
	mdouble endcap_delta_cone[2];	//< endcap inner radius cone offset
	mdouble endcap_mu[2];			//< permeability for endcap
	unsigned int eSegs[2];			//< endcap segments
};

/// Data-sampling cell
class fieldCell {
public:
	/// constructor
	fieldCell(vec3 l, vec3 u, unsigned int n1, unsigned int n2, unsigned int n3): ll(l), ur(u), nx(n1), ny(n2), nz(n3), vx(5), vy(5), vz(5) {}
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
