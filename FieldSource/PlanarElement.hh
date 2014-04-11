/* 
 * PlanarElement.hh, part of the RotationShield program
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

#ifndef PLANARELEMENT_HH
/// Make sure this header is only loaded once
#define PLANARELEMENT_HH

#include "ReactiveElement.hh"

/// a ReactiveElement on a plane
class PlanarElement: public ReactiveElement {
public:
	/// constructor
	PlanarElement(Plane pl): ReactiveElement(), p(pl) { }
		
	/// replicate around a new plane
	virtual PlanarElement* replicate(Plane pl) const = 0;
	
	/// replicate with new annulus specification
	virtual PlanarElement* reference(annulusSpec a) const { return replicate(Plane(a)); }
	
	/// replicate rotated around z axis
	virtual PlanarElement* replicateRotated(double th) const = 0;
	
	/// Visualize the element
	virtual void _visualize() const;
	
	Plane p;		///< The plane in which the element resides
};


#endif
