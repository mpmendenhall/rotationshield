#ifndef BOXEL_HH
#define BOXEL_HH 1

#include "CompoundPlane.hh"

class Boxel: public CompoundPlane {
public:
	/// constructor
	Boxel(annulusSpec a, mdouble mu, mdouble thickness);
	
	/// destructor
	virtual ~Boxel() {}
	
	/// generate initial reference element
	virtual PlanarElement* reference(annulusSpec a) const { Boxel* w = new Boxel(a,murel,d); w->preSolve(); return w; }
	
	/*
	/// Visualize the element
	virtual void visualize(bool top = true, mdouble logmax = 3.0) const {
		if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); }
		for(unsigned int i=0; i<drawSurfaces; i++)
			els[i]->visualize(false);
		if(top) vsr::Visr::W->stopRecording();
	}
	*/
	
	/// replicate reference element to other angle
	virtual PlanarElement* replicateRotated(mdouble th) const {
		Boxel* C = new Boxel(p.zrotated(th), elResponse);
		for(unsigned int i=0; i<els.size(); i++)
			C->append(els[i]->replicateRotated(th));
		return C;
	}
	
protected:
	
	/// private constructor for replication
	Boxel(Plane pl, const std::vector<mmat>& R): CompoundPlane(pl,R) {}
	
	mdouble murel;	//< relative permeability
	mdouble d;		//< thickness		
	static unsigned int drawSurfaces; //< number of surfaces to draw
};

#endif
