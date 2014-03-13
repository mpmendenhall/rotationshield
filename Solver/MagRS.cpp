#include "MagRS.hh"

BField_Protocol* BField_Protocol::BFP = new BField_Protocol();

void MagRSCombiner::addSet(ReactiveSet* R) {
	ReactiveSetCombiner::addSet(R);
};

void MagRSCombiner::calculateIncident(const FieldSource& f) {
	for(std::vector<ReactiveSet*>::iterator it = mySets.begin(); it != mySets.end(); it++) {
		MagF_Responder* MR = dynamic_cast<MagF_Responder*>(*it);
		MR->calculateIncident(f);
	}
}

vec3 MagRSCombiner::fieldAt(const vec3& v) const {
	vec3 B;
	for(std::vector<ReactiveSet*>::const_iterator it = mySets.begin(); it != mySets.end(); it++) {
		FieldSource* FS = dynamic_cast<FieldSource*>(*it);
		B += FS->fieldAt(v);
	}
	return B;
}

void MagRSCombiner::_visualize() const {
	for(std::vector<ReactiveSet*>::const_iterator it = mySets.begin(); it != mySets.end(); it++) {
		FieldSource* FS = dynamic_cast<FieldSource*>(*it);
		FS->_visualize();
	}
}
