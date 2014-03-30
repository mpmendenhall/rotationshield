/* 
 * RefCounter.cpp, part of the RotationShield program
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

#include "RefCounter.hh"
#include <stdlib.h>

int RefCounter::nReferencedObjects = 0;
int RefCounter::nTotalReferences = 0;

void RefCounter::release() const {
	if(!refcount) {
		printf("Object '%s' released with reference count of 0! I die!\n",ref_name.c_str());
		exit(-1);
	}
	--refcount;
	--nTotalReferences;
	if(!refcount) delete(this);
}

RefCounter::~RefCounter() {
	if(refcount) {
		printf("Object '%s' deleted with reference count of %i! I die!\n",ref_name.c_str(),refcount);
		exit(-1);
	}
	--nReferencedObjects;
}
