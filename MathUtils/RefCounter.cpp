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