#ifndef GRIDBOXEL_HH
#define GRIDBOXEL_HH 1

#include "CompoundElement.hh"
#include "Planesource.hh"

class GridBoxel: public CompoundElement {
public:
	/// constructor
	GridBoxel(annulusSpec a, mdouble mu, mdouble thickness, unsigned int nz, unsigned int nt, unsigned int nd);
	/// constructor for replication
	
	/// destructor
	virtual ~GridBoxel() {}
	
	/// generate initial reference element
	virtual ReactiveElement* reference(annulusSpec a) const { GridBoxel* w = new GridBoxel(a,murel,d,nZ,nT,nD); w->preSolve(); return w; }
	
	/// Visualize the element
	virtual void _visualize() const {
		unsigned int elmax = els.size();
		if(drawSurfaces<3)
			elmax = 2*nZ*nT;
		for(unsigned int i=0; i<elmax; i++)
			els[i]->_visualize();
	}
	
	/// replicate reference element to other angle
	virtual ReactiveElement* replicateRotated(mdouble th) const {
		GridBoxel* C = new GridBoxel(p.zrotated(th), elResponse, nZ, nT, nD);
		for(unsigned int i=0; i<els.size(); i++)
			C->append(els[i]->replicateRotated(th));
		return C;
	}
	
protected:
	
	/// private constructor for replication
	GridBoxel(Plane pl, const std::vector<mat3>& R, unsigned int nz, unsigned int nt, unsigned int nd): CompoundElement(pl,R), nZ(nz), nT(nt), nD(nd) {}
	
	mdouble murel;	//< relative permeability
	mdouble d;		//< thickness		
	unsigned int nZ;
	unsigned int nT;
	unsigned int nD;
	static unsigned int drawSurfaces; //< number of surfaces to draw
};

#endif
