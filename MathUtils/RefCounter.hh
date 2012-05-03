#ifndef REFCOUNTER_HH
#define REFCOUNTER_HH 1

#include <cassert>
#include <stdio.h>

class RefCounter {
public:
	RefCounter(): refcount(0) { ++nReferencedObjects; }
	virtual ~RefCounter() { assert(refcount==0); --nReferencedObjects; }
	void retain() const { ++refcount; ++nTotalReferences; }
	void release() const { assert(refcount>0); --refcount; --nTotalReferences; if(!refcount) delete(this); }
	static int referencedObjectCount() { return nReferencedObjects; }
	static int totalReferenceCount() { return nTotalReferences; }
protected:
	mutable unsigned int refcount;
	static int nReferencedObjects;
	static int nTotalReferences;
};


class VerboseRefCounter {
public:
	VerboseRefCounter(): refcount(0) { printf("*%p created\n",this); }
	virtual ~VerboseRefCounter() { assert(refcount==0);  printf("*%p destroyed\n",this); }
	void retain() const { ++refcount;  printf("*%p retained [%i]\n",this,refcount); }
	void release() const { assert(refcount>0); --refcount;  printf("*%p released [%i]\n",this,refcount); if(!refcount) delete(this); }
protected:
	mutable unsigned int refcount;
};


#endif
