#ifndef REACTIVEELEMENT_HH
#define REACTIVEELEMENT_HH 1

#include "FieldSource.hh"
#include "Color.hh"
#include <cmath>
#include "VarMat.hh"

/// FieldSource that responds linearly to incident fields in order to satisfy boundary conditions
class ReactiveElement: public FieldSource {
public:
	
	/// constructor
	ReactiveElement(): state(mvec()) { }

	/// number of degrees of freedom
	virtual unsigned int nDF() const = 0;
	
	/// destructor
	virtual ~ReactiveElement() {}
		
	/// get fields at a point due to state components
	virtual mmat fieldAtComponents(vec3 p0) const = 0;
	/// resulting state in response to incident field f
	virtual mvec responseToFieldSource(const FieldSource* f) const = 0;
	/// interaction matrix with another ReactiveElement
	virtual mmat interactionWith(const ReactiveElement* e) const = 0;
	/// set state in reaction to applied field
	void reactTo(const FieldSource* f) { setState(responseToFieldSource(f)); }
	
	/// get field produced at a point
	vec3 fieldAt(const vec3& v) const { 
		mvec f =  fieldAtComponents(v)*state;
		return vec3(f[0],f[1],f[2]);
	}

	/// get state vector
	const mvec& getState() const { return state; }
	
	/// set state vector
	virtual void setState(mvec v) { assert(v.size() == nDF()); state = v; }
		
protected:
	
	mvec state;		//< vector describing state of element (e.g. surface current for a PlaneSource)
};

#endif
