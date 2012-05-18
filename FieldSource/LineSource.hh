#ifndef LINESOURCE_HH
#define LINESOURCE_HH 1

#include "FieldSource.hh"


/// Magenetic field source from a line of current (i.e. a straight segment of wire)
class LineSource: public FieldSource {
public:
	/// Constructor
	LineSource(vec3 startv, vec3 endv, mdouble current): FieldSource(), l(startv,endv), j(current) {}
	/// Destructor
	~LineSource() {}
	
	/// Field at a specified point
	vec3 fieldAt(const vec3& v) const;
	/// Averages the field from a LineSource over a given Line, with special cases for parallel and perpendicular lines
	virtual vec3 fieldOverLine(Line l) const;
	
	/// Get the LineSource's Line #l
	Line getLine() { return l; }
	/// Total vector current of the LineSource
	vec3 dipole() const { return (l.e- l.s)*j; }
	/// Print info to stdout
	void display() const { printf("Linesource (j=%g):\n\t",(double)j); l.display(); }
	/// Visualize the field source
	virtual void visualize(bool top = true, mdouble scaling = 1.0) const {
		if(top) { vsr::startRecording(true); vsr::clearWindow(); } 
		//l.visualize(false); 
		l.visualizeDirected(sign(j));
		if(top) vsr::stopRecording(); 
	}
	
	/// Integrates \f$ \int_0^x \frac{(t^2+a^2)^{1/2}}{(t^2+b^2)^{3/2}}dt \f$ (\f$ a>b \f$)
	/** See Abramowitz and Stegun 17.4.41 */
	static mdouble sF1(mdouble a2, mdouble b2, mdouble x);
	/// Integrates \f$ \int_0^x \frac{1}{(t^2+a^2)^{3/2}(t^2+b^2)^{1/2}}dt \f$
	/** See Abramowitz and Stegun 17.4.51 */
	static mdouble sF2(mdouble a2, mdouble b2, mdouble x);
	/// Integrates \f$ \int_{x_1}^{x_2} \frac{(t^2+a^2)^{1/2}}{(t^2+b^2)^{1/2}}dt \f$ (\f$ a>b \f$)
	static mdouble sF1(mdouble a2, mdouble b2, mdouble x1, mdouble x2) { return sF1(a2,b2,x2) - sF1(a2,b2,x1); }
	/// Integrates \f$ \int_{x_1}^{x_2} \frac{1}{(t^2+a^2)^{3/2}(t^2+b^2)^{1/2}}dt \f$
	static mdouble sF2(mdouble a2, mdouble b2, mdouble x1, mdouble x2) { return sF2(a2,b2,x2) - sF2(a2,b2,x1); }
	
private:
	const Line l;	//< The line along which current flows
	const mdouble j; //< The current
};

#endif
