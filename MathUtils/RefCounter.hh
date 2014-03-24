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
