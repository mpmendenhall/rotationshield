/* 
 * MagRS.cpp, part of the RotationShield program
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

#include "MagRS.hh"

BField_Protocol* BField_Protocol::BFP = new BField_Protocol();

void MagRSCombiner::addSet(ReactiveSet* R) {
	ReactiveSetCombiner::addSet(R);
}

void MagRSCombiner::calculateIncident(const FieldSource& f) {
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++) {
		MagF_Responder* MR = dynamic_cast<MagF_Responder*>(*it);
		MR->calculateIncident(f);
	}
}

vec3 MagRSCombiner::fieldAt(const vec3& v) const {
	vec3 B;
	for(std::vector<ReactiveSet*>::const_iterator it = mySets.begin(); it != mySets.end(); it++) {
		FieldSource* FS = dynamic_cast<FieldSource*>(*it);
		B += FS->fieldAt(v);
	}
	return B;
}

void MagRSCombiner::_visualize() const {
	for(std::vector<ReactiveSet*>::const_iterator it = mySets.begin(); it != mySets.end(); it++) {
		FieldSource* FS = dynamic_cast<FieldSource*>(*it);
		FS->_visualize();
	}
}
