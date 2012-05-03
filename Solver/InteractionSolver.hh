#ifndef INTERACTIONSOLVER_HH
#define INTERACTIONSOLVER_HH 1

#include "ReactiveElement.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "MiscUtils.hh"
#include <stdlib.h>
#include <math.h>

/// virtual base class for solving interacting systems
class InteractionSolver: public FieldSource {
public:
	/// Constructor
	InteractionSolver(unsigned int gs = 1);
	/// Destructor
	virtual ~InteractionSolver();
	/// add an interacting element
	void addSurfacel(ReactiveElement* e);
	/// get number of interacting elements
	unsigned int N() const { return surfacels.size(); }
	/// Solve for the Greene's Function of the system
	virtual void solve(bool checkinversion = false) = 0;
	
	/// Calculate the #incidentState on the shield produced by the specified FieldSource
	void calculateIncident(FieldSource* f);	
	/// Magnetic field at v produced by the shield elements
	vec3 fieldAt(const vec3& v) const;	
	/// Magnetic field at v produced by a subset of the field elements
	vec3 partialField(const vec3& v, unsigned int n1, int nf = -1) const;	
	/// Calculate the shield's #finalState in response to the #incidentState
	virtual void calculateResult() = 0;
	
	/// visualization of the interacting elements
	virtual void visualize(bool top = true, mdouble scaling = 1.0) const;	
	/// get final state after calculations
	mvec getFinalState(unsigned int i) const;	
		
	mvec incidentState; 	//< surfacels' state in response to incident field
	mvec finalState;		//< surfacels' final state after interactions
	
protected:
	
	/// Interaction between two shield elements
	mmat getInteraction(int x0, int x1) const;
	/// sets surfacels to final state
	void setSurfacelsOutputState() { for(unsigned int i=0; i<N(); i++) surfacels[i]->setState(getFinalState(i)); }
	/// get array index for df of given element number
	unsigned int index(unsigned int el, unsigned int df) const { return groupIndex[el/groupSize]+df*groupSize+(el%groupSize); }
	
	unsigned int ndf;							//< total number of degrees of freedom
	const unsigned int groupSize;				//< elements will be loaded in groups of groupSize elements
	bool groupsComplete;						//< whether last group is completely filled
	FieldSource* fInc;							//< incident field source
	std::vector<ReactiveElement*> surfacels;	//< the shield surface elements
	std::vector<unsigned int> groupIndex;		//< starting index for each group
	bool verbose;								//< whether to print progress status
};

#endif
