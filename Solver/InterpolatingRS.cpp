#include "InterpolatingRS.hh"
#include <cassert>

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
	unsigned int p,z,d;
	DF_address(DF,p,z,d);
	assert(d<G.size() && p<nPhi && z<nZ);
	G[d]->set(z,p,v);
}
