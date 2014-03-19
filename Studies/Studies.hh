/// \file Studies.hh \brief Utility classes for setting up Cos Theta Coil studies

#ifndef STUDIES_HH
/// Makes sure to only load this file once
#define STUDIES_HH 1

#include "CosThetaBuilder.hh"
#include "MixedSource.hh"
#include "SymmetricSolver.hh"
#include "SurfaceProfiles.hh"
#include "MagRS.hh"
#include "QFile.hh"
#include <vector>

/// geometry reference points for shield segments
enum GeomRefPt {
	GEOMREF_ORIGIN,		//< absolute origin
	GEOMREF_CENTER,		//< center of frame
	GEOMREF_NEGAXIS,	//< negative z, r = 0
	GEOMREF_NEGRADIUS,	//< negative z, r = r0
	GEOMREF_POSRADIUS,	//< positive z, r = r0
	GEOMREF_POSAXIS		//< positive z, r = r0
};

/// section of a shield
class shieldSection {
public:
	/// constructor
	shieldSection();
	
	/// shield info as Stringmap
	Stringmap getInfo() const;
	
	mdouble mu;				//< material relative permeability
	unsigned int cSegs;		//< constant-width z segments
	unsigned int vSegs;		//< variable z segments
	GeomRefPt endpts[2];	//< endpoint reference locations
	vec2 endoff[2];			//< endpoint offsets from reference locations
};

/// Shield surface framework, defining reference points around which surface segments are built
class shieldFrame {
public:
	/// constructor
	shieldFrame(): pSegs(0), length(0), radius(0), RSC(NULL) {}
	/// destructor
	virtual ~shieldFrame() { if(RSC) delete RSC; }
	
	unsigned int pSegs;						//< phi segments
	mdouble length;							//< shield length
	mdouble radius;							//< shield radius
	std::vector<shieldSection> mySections;	//< sections of shield
	MagRSCombiner* RSC;						//< shield system to be solved
	SymmetricSolver SS;						//< solver for shield
	
	/// return reference point location
	vec2 refPt(GeomRefPt p) const;
	/// construct into mixed source
	void construct(MixedSource& ms, CosThetaBuilder* ct, const std::string& fcache = "");
	/// shield info as Stringmap
	Stringmap getInfo() const;
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
	
	CosThetaBuilder coil;	//< coil providing measured field
	shieldFrame shield;		//< shielding around coil
	fieldCell cell;			//< measurement cell
	
	bool saveGrid;			//< whether to save the grid data to file
	std::string basedir;	//< directory for output
	
protected:
	MixedSource* ms;		//< combined coil+shield field source
	
};

#endif
