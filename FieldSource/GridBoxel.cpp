#include "GridBoxel.hh"
#include <vector>

unsigned int GridBoxel::drawSurfaces = 6;

GridBoxel::GridBoxel(annulusSpec a, mdouble mu, mdouble thickness, unsigned int nz, unsigned int nt, unsigned int nd):
CompoundElement(a), murel(mu), d(thickness), nZ(nz), nT(nt), nD(nd) {
	
	std::vector<annulusSpec>::iterator it;
	std::vector<Plane>::iterator itp;
	std::vector<annulusSpec> slices = a.subdivide(nZ,nT);
	
	std::vector<AnnularBoxel> chunks = std::vector<AnnularBoxel>();
	for(it = slices.begin(); it != slices.end(); it++) {
		chunks.push_back(AnnularBoxel(*it,d));
		append(new PlaneSource(Plane(chunks.back().top),murel));
		append(new PlaneSource(Plane(chunks.back().bottom),murel));
	}
	
	for(unsigned int t=0; t<nT; t++) {
		std::vector<annulusSpec> v1 = chunks[t].front.subdivide(nD,1);
		for(it = v1.begin(); it != v1.end(); it++)
			append(new PlaneSource(Plane(*it),murel));
		std::vector<annulusSpec> v2 = chunks[nZ*nT-1-t].back.subdivide(nD,1);
		for(it = v2.begin(); it != v2.end(); it++)
			append(new PlaneSource(Plane(*it),murel));
	}
	
	for(unsigned int z=0; z<nZ; z++) {
		std::vector<Plane> v1 = chunks[z*nT].left.subdivide(nD,1);
		for(itp = v1.begin(); itp != v1.end(); itp++)
			append(new PlaneSource(*itp,murel));
		std::vector<Plane> v2 = chunks[z*nT+nT-1].right.subdivide(nD,1);
		for(itp = v2.begin(); itp != v2.end(); itp++)
			append(new PlaneSource(*itp,murel));
	}
}