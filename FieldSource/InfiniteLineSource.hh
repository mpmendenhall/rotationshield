/* 
 * InfiniteLineSource.hh, part of the RotationShield program
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

#ifndef INFINITELINESOURCE_HH
/// Make sure this header is only loaded once
#define INFINITELINESOURCE_HH

#include "FieldSource.hh"


/// Magenetic field source from a line of current (i.e. a straight segment of wire)
class InfiniteLineSource: public FieldSource {
public:
	/// Constructor
	InfiniteLineSource(vec3 startv, vec3 endv, double current): FieldSource(), l(startv,endv), j(current) {}
	/// Constructor for 2D positioning
	InfiniteLineSource(vec2 v0, double current): FieldSource(), l(vec3(v0[0],v0[1],-0.5),vec3(v0[0],v0[1],0.5)), j(current) {}
	/// Destructor
	~InfiniteLineSource() {}
	
	/// Field at a specified point
	vec3 fieldAt(const vec3& v) const;
	/// Averages the field from a LineSource over a given Line, with special cases for parallel and perpendicular lines
	//virtual vec3 fieldOverLine(Line l) const;
	
	/// Print info to stdout
	void display() const { printf("InfiniteLinesource (j=%g):\n\t",(double)j); l.display(); }
	/// Visualize the field source
	virtual void _visualize() const { l.visualizeDirected(sign(j)); }
		
private:
	const Line l;	///< The line along which current flows
	const double j; ///< The current
};

#endif
