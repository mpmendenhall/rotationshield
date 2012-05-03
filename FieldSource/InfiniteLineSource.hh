#ifndef INFINITELINESOURCE_HH
#define INFINITELINESOURCE_HH 1

#include "FieldSource.hh"


/// Magenetic field source from a line of current (i.e. a straight segment of wire)
class InfiniteLineSource: public FieldSource {
public:
	/// Constructor
	InfiniteLineSource(vec3 startv, vec3 endv, mdouble current): FieldSource(), l(startv,endv), j(current) {}
	/// Constructor for 2D positioning
	InfiniteLineSource(vec2 v0, mdouble current): FieldSource(), l(vec3(v0[0],v0[1],-0.5),vec3(v0[0],v0[1],0.5)), j(current) {}
	/// Destructor
	~InfiniteLineSource() {}
	
	/// Field at a specified point
	vec3 fieldAt(const vec3& v) const;
	/// Averages the field from a LineSource over a given Line, with special cases for parallel and perpendicular lines
	//virtual vec3 fieldOverLine(Line l) const;
	
	/// Print info to stdout
	void display() const { printf("InfiniteLinesource (j=%g):\n\t",(double)j); l.display(); }
	/// Visualize the field source
	virtual void visualize(bool top = true, mdouble scaling = 1.0) const {
		if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); } 
		//l.visualize(false); 
		l.visualizeDirected(sign(j));
		if(top) vsr::Visr::W->stopRecording(); 
	}
		
private:
	const Line l;	//< The line along which current flows
	const mdouble j; //< The current
};

#endif
