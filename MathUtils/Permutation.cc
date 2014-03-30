/* 
 * Permutation.cc, part of the RotationShield program
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

#include "Permutation.hh"

Permutation::Permutation(unsigned int n): data(std::vector<unsigned int>(n)) {
	for(unsigned int i=0; i<n; i++)
		data[i] = i;
}

Permutation& Permutation::nshuffle(int n) {
	for(unsigned int i=0; i<size()/n; i++) 
		for(int j=0; j<n; j++)
			data[j*size()/n+i] = i*n+j;
	return *this;
}

void Permutation::swap(unsigned int a, unsigned int b) {
	std::swap(data[a],data[b]);
}

Permutation Permutation::inverse() const {
	Permutation inv = Permutation(size());
	for(unsigned int i=0; i<size(); i++) inv[data[i]] = i;
	return inv;
}

Permutation& Permutation::invert() {
	std::vector<unsigned int> newdat = std::vector<unsigned int>(size());
	for(unsigned int i=0; i<size(); i++) newdat[data[i]] = i;
	data = newdat;
	return *this;
}

const Permutation Permutation::operator*(const Permutation& p) const {
	assert(p.size() == size());
	Permutation o = Permutation(size());
	for(unsigned int i=0; i<size(); i++) o[i] = p[data[i]];
	return o;
}

int main() {
	return 0;
}