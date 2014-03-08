#ifndef MIXEDSOURCE_HH
#define MIXEDSOURCE_HH 1

#include "FieldSource.hh"
#include "LineSource.hh"
#include <vector>

/// Combines several FieldSource objects into one (e.g. many LineSource wires into a Cos Theta coil)
class MixedSource: public FieldSource {
public:
	/// Constructor
	MixedSource(): FieldSource(), displayColor(vec3(1.0,0.0,0.0)), sources(std::vector<const FieldSource*>()) {}
	
	/// Destructor
	virtual ~MixedSource() {
		for(unsigned int i=0; i<sources.size(); i++)
			sources[i]->release();
	}
	
	/// Field at a specified point
	vec3 fieldAt(const vec3& v) const;
	/// Field averaged over the given Line
	virtual vec3 fieldOverLine(Line l) const;
	/// Field averaged over the given Plane
	virtual vec3 fieldOverPlane(Plane p) const;
	/// Add a new FieldSource to the source currents
	void addsource(const FieldSource* fs) { fs->retain(); sources.push_back(fs); }
	/// Add FieldSource's specified in a file
	void loadSourcesFile(FILE* f, mdouble scale);
	
	/// Net current of the arrangement
	vec3 dipole() const { vec3 d = vec3(); for(unsigned int i=0; i<sources.size(); i++) d += sources[i]->dipole(); return d; }
	
	/// Add an arc of current (approximated by many LineSource segments)
	void arc(vec3 start, vec3 end, mdouble j, int nsegs = 1);
	/// Add a loop of current (approximated by many LineSource segments)
	void loop(mdouble z, mdouble r, mdouble j, int nsides = 32);
	
	/// Print info to stdout
	void display() const { printf("Multisource:\n"); for(unsigned int i=0; i<sources.size(); i++) sources[i]->display(); }
	/// Visualize the field source
	virtual void visualize(bool top = true, mdouble scaling = 1.0) const {
		if(top) { vsr::startRecording(true); vsr::clearWindow(); }
		vsr::setColor(displayColor[0],displayColor[1],displayColor[2]);
		for(unsigned int i=0; i<sources.size(); i++) sources[i]->visualize(false, scaling);
		if(top) vsr::stopRecording();
	}
	
	vec3 displayColor;
	
private:

	std::vector<const FieldSource*> sources; //< component sources
};

#endif
