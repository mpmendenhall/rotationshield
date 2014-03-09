#ifndef COMPOUNDELEMENT_HH
#define COMPOUNDELEMENT_HH 1

#include "PlanarElement.hh"
#include "GenericSolver.hh"
#include "UniformField.hh"
#include <vector>

/// Collection of several ReactiveElements that act as one unit
class CompoundPlane: public PlanarElement {
public:	
	/// constructor
	CompoundPlane(Plane pl):
	PlanarElement(pl), els(std::vector<PlanarElement*>()), elResponse(std::vector<mmat>()) { rmat = p.projectionMatrix(); }
	
	/// destructor
	virtual ~CompoundPlane() { while(els.size()) { els.back()->release(); els.pop_back(); } }
	
	/// interaction matrix with another ReactiveElement
	virtual mmat interactionWithComponents(const ReactiveElement* e) const { return rmat * e->fieldAtComponents(p.o); }
	
	/// resulting state in response to incident field f
	virtual mvec responseToFieldSource(const FieldSource* f) const { return rmat * f->fieldOverPlane(p); }
	
	/// append a new subelement
	void append(PlanarElement* e) { e->retain(); els.push_back(e); }
	
	/// generate initial reference element
	virtual PlanarElement* reference(annulusSpec a) const { assert(false); return NULL; }
	
	/// return number of external DF (always 3)
	virtual unsigned int nDF() const { return 3; }
	
	/// replicate reference element to other angle
	virtual PlanarElement* replicateRotated(mdouble th) const {
		CompoundPlane* C = new CompoundPlane(p.zrotated(th), elResponse);
		for(unsigned int i=0; i<els.size(); i++)
			C->append(els[i]->replicateRotated(th));
		return C;
	}
	
	/// replicate around a new plane
	virtual PlanarElement* replicate(Plane pl) const { assert(false); return NULL; }
	
	/// get fields due to state components
	virtual mmat fieldAtComponents(vec3 p0) const {
		mmat M = mmat(3,nDF());
		for(unsigned int i=0; i<els.size(); i++)
			M += els[i]->fieldAtComponents(p0) *  elResponse[i];
		return M;
	}
	
	/// get field at given point
	virtual vec3 fieldAt(const vec3& v) const {
		vec3 f = vec3(); 
		for(unsigned int i=0; i<els.size(); i++)
			f += els[i]->fieldAt(v);
		return f;
	}

	
	/// set element state => subelement states
	virtual void setState(mvec v) {
		assert(v.size() == nDF());
		state = v;
		for(unsigned int i=0; i<els.size(); i++)
			els[i]->setState(elResponse[i]*v);
	}
	
	/*
	/// Visualize the element's subelements
	virtual void visualize(bool top = true, mdouble logmax = 3.0) const {
		if(top) { vsr::startRecording(); vsr::clearWindow(); }
		for(unsigned int i=0; i<els.size(); i++)
			els[i]->visualize(false);
		if(top) vsr::stopRecording();
	}
	*/
	
	/// solve for elements' response to external fields
	void preSolve() {
		
		printf("Pre-solving for subassembly (%i components)...\n",(int)els.size());
		
		GenericSolver* GS = new GenericSolver();
		for(unsigned int i=0; i<els.size(); i++)
			GS->addSurfacel(els[i]);
		//addAuxElements(GS);
		GS->solve();
		
		elResponse = std::vector<mmat>(els.size());
		for(unsigned int i=0; i<els.size(); i++)
			elResponse[i] = mmat(els[i]->nDF(),nDF());
		
		// calculate response to various initial fields
		for(unsigned int axis = 0; axis<nDF(); axis++) {
			vec3 incident = vec3();
			for(unsigned int i=0; i<3; i++)
				incident[i] = rmat(axis,i);
			GS->calculateIncident(new UniformField(incident));
			GS->calculateResult();
			for(unsigned int i=0; i<els.size(); i++) {
				mvec f = GS->getFinalState(i);
				for(unsigned int j=0; j<els[i]->nDF(); j++)
					elResponse[i](j,axis) = f[j];
			}
		}
	}
	
	
	
protected:	
	
	/// private constructor for replication
	CompoundPlane(Plane pl, const std::vector<mmat>& R):
	PlanarElement(pl), els(std::vector<PlanarElement*>()), elResponse(R) { rmat = p.projectionMatrix(); }
	
	/// add extra elements if needed when calculating response
	//virtual void addAuxElements(GenericSolver* GS) const {}
	
	std::vector<PlanarElement*> els;	//< constituent subelements
	std::vector<mmat> elResponse;		//< elements' response to applied fields in p frame
	mmat rmat;							//< rotation matrix into local frame
};

#endif
