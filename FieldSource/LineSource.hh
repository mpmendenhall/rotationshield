/* 
 * LineSource.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef LINESOURCE_HH
#define LINESOURCE_HH 1

#include "FieldSource.hh"


/// Magenetic field source from a line of current (i.e. a straight segment of wire)
class LineSource: public FieldSource {
public:
	/// Constructor
	LineSource(vec3 startv, vec3 endv, double current): FieldSource(), l(startv,endv), j(current) {}
	/// Destructor
	~LineSource() {}
	
	/// Field at a specified point
	vec3 fieldAt(const vec3& v) const;
	/// Averages the field from a LineSource over a given Line, with special cases for parallel and perpendicular lines
	virtual vec3 fieldOverLine(Line l) const;
	
	/// Get the LineSource's Line #l
	Line getLine() { return l; }
	/// Total vector current of the LineSource
	vec3 net_current() const { return (l.e- l.s)*j; }
	/// Print info to stdout
	void display() const { printf("Linesource (j=%g):\n\t",(double)j); l.display(); }
	/// Visualize the field source
	virtual void _visualize() const { l.visualizeDirected(sign(j)); }
	
	/// Integrates \f$ \int_0^x \frac{(t^2+a^2)^{1/2}}{(t^2+b^2)^{3/2}}dt \f$ (\f$ a>b \f$)
	/** See Abramowitz and Stegun 17.4.41 */
	static double sF1(double a2, double b2, double x);
	/// Integrates \f$ \int_0^x \frac{1}{(t^2+a^2)^{3/2}(t^2+b^2)^{1/2}}dt \f$
	/** See Abramowitz and Stegun 17.4.51 */
	static double sF2(double a2, double b2, double x);
	/// Integrates \f$ \int_{x_1}^{x_2} \frac{(t^2+a^2)^{1/2}}{(t^2+b^2)^{1/2}}dt \f$ (\f$ a>b \f$)
	static double sF1(double a2, double b2, double x1, double x2) { return sF1(a2,b2,x2) - sF1(a2,b2,x1); }
	/// Integrates \f$ \int_{x_1}^{x_2} \frac{1}{(t^2+a^2)^{3/2}(t^2+b^2)^{1/2}}dt \f$
	static double sF2(double a2, double b2, double x1, double x2) { return sF2(a2,b2,x2) - sF2(a2,b2,x1); }
	
private:
	const Line l;	//< The line along which current flows
	const double j; //< The current
};

#endif
