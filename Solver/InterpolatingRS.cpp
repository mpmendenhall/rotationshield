#include "InterpolatingRS.hh"
#include <cassert>
#include <stdio.h>

void InterpolatingRS2D::clear_data() {
	for(unsigned int i=0; i<G.size(); i++)
		delete G[i];
	G.clear();
}

void InterpolatingRS2D::make_grids(unsigned int nz, unsigned int ndf) {
	clear_data();
	nDFi = ndf;
	nZ = nz;
	for(unsigned int i=0; i<nDFi; i++) {
		G.push_back(new BicubicGrid(nZ,nPhi));
		if(nPhi>1) {
			G.back()->bc[1] = IB_CYCLIC;
			G.back()->setUserRange(0,1,true,0.5);
			G.back()->setUserRange(0,1,false,0.5);
		}
	}
}

void InterpolatingRS2D::DF_address(unsigned int DF, unsigned int& p, unsigned int& z, unsigned int& d) const {
	p = DF%nPhi;
	z = (DF/nPhi)%nZ;
	d = DF/(nPhi*nZ);
}

void InterpolatingRS2D::_setDF(unsigned int DF, double v) {
	assert(DF < nDF());
	assert(v==v && fabs(v)<1e6);
	unsigned int p,z,d;
	DF_address(DF,p,z,d);
	assert(d<G.size() && p<nPhi && z<nZ);
	G[d]->set(z,p,v);
}

mvec InterpolatingRS2D::interpl_DF(vec2 l) const {
	mvec v(nDFi);
	for(unsigned int i=0; i<nDFi; i++)
		v[i] = (*G[i])(l[0],l[1]);
	return v;
}

void InterpolatingRS2D::setToroidal() {
	isToroidal = true;
	for(unsigned int i=0; i<nDFi; i++)
		G[i]->bc[0] = IB_CYCLIC;
}

void InterpolatingRS2D::printData() const {
	for(unsigned int i=0; i<nDFi; i++) {
		printf("Data table for layer %i:\n",i);
		G[i]->printData();
	}
}
