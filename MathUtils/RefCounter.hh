/* 
 * RefCounter.hh, part of the RotationShield program
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

#ifndef REFCOUNTER_HH
#define REFCOUNTER_HH 1

#include <cassert>
#include <stdio.h>
#include <string>

class RefCounter {
public:
	RefCounter(const std::string& nm = "RefCounter"): ref_name(nm), refcount(0) { ++nReferencedObjects; }
	virtual ~RefCounter();
	virtual void retain() const { ++refcount; ++nTotalReferences; }
	virtual void release() const;
	static int referencedObjectCount() { return nReferencedObjects; }
	static int totalReferenceCount() { return nTotalReferences; }
	
	std::string ref_name;	//< name for this object, for tracing down deletion issues
	
protected:
	mutable unsigned int refcount;
	static int nReferencedObjects;
	static int nTotalReferences;
};


class VerboseRefCounter: public RefCounter {
public:
	VerboseRefCounter(): RefCounter() { printf("'%s' *%p created\n",ref_name.c_str(),(void*)this); }
	virtual ~VerboseRefCounter() { printf("'%s' *%p destroyed\n",ref_name.c_str(),(void*)this); }
	virtual void retain() const { RefCounter::retain();  printf("'%s' *%p retained [%i]\n",ref_name.c_str(),(void*)this,refcount); }
	virtual void release() const { printf("'%s' *%p released [%i]\n",ref_name.c_str(),(void*)this,refcount-1); RefCounter::release(); }
};


#endif
