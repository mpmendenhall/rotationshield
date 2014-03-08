/// \file ReactiveSet.hh \brief Virtual base class for collections of interacting degrees of freedom

#ifndef REACTIVESET_HH
/// Makes sure to only load this file once
#define REACTIVESET_HH 1

#include "VarMat.hh"
#include "VarVec.hh"
#include <vector>
#include <cassert>

typedef VarMat<mdouble> mmat;
typedef VarVec<mdouble> mvec;

/// a set of linearly interacting degrees of freedom
class ReactiveSet {
public:
	/// constructor
	ReactiveSet(unsigned int nph=1): nPhi(nph), myClass(RS_BASE) {}
	/// destructor
	virtual ~ReactiveSet() {}
	
	/// total number of degrees of freedom
	virtual unsigned int nDF() const = 0;
	
	/// reset interaction term counter
	virtual void startInteractionScan() = 0;
	/// possibly arbitrarily ordered interaction matrix entries
	virtual mdouble nextInteractionTerm(unsigned int& i, unsigned int& j) = 0;
	
	const unsigned int nPhi;	//< internal periodic symmetry

	enum RS_CLASSTYPE {
		RS_BASE = 1<<0,
		RS_FS = 1<<1,
		RS_SS = 1<<2
	} myClass;	//< class indicator for specialized interactions
	
	mvec incidentState; 	//< non-interacting initial state vector
	/// set (final) state
	virtual void setFinalState(const mvec& v) { finalState = v; }

protected:

	mvec finalState;					//< final state after interactions
	mutable unsigned int rterm;			//< current interaction term
	
};

/// linearly interacting DF grouped into logical sub-units
class ReactiveUnitSet: public ReactiveSet {
public:
	/// constructor
	ReactiveUnitSet(unsigned int nph=1): ReactiveSet(nph) { group_start.push_back(0); }
	
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return group_start.back(); }
	
	/// reset interaction term counter
	virtual void startInteractionScan() { ic_i = ic_j = ic_di = ic_dj = 0; ic_m = interactionBetween(ic_i, ic_j); }
	/// possibly arbitrarily ordered interaction matrix entries
	virtual mdouble nextInteractionTerm(unsigned int& i, unsigned int& j);
	
	
	/// set state for i^th sub-element
	virtual void setState(unsigned int i, const mvec& v) { assert(false); }
	/// extracts individual element final state from state vector
	mvec getFinalState(unsigned int i) const { assert(false); }
	
protected:
	
	/// get DF index for sub-element DF
	unsigned int df_subindex(unsigned int el, unsigned int df) const { return group_start[el/nPhi] + df*nPhi + (el%nPhi); }
	
	/// add a group of sub-elements with given number of DF
	void add_DF_group(unsigned int N) { group_DF.push_back(N); group_start.push_back(group_start.back()+nPhi*N); }
	
	/// interaction matrix between subunits i and j
	virtual mmat interactionBetween(unsigned int i, unsigned int j) const = 0;
	
	// interaction calculator cache variables
	mmat ic_m;			//< cached interaction matrix between sub-elements
	unsigned int ic_i;	//< i element interaction
	unsigned int ic_j;	//< j element interaction
	unsigned int ic_di;	//< i element DF
	unsigned int ic_dj;	//< j element DF
	
	unsigned int n_subels() const { return group_DF.size()*nPhi; }
	std::vector<unsigned int> group_DF;		//< number of DF for elements in this group
	std::vector<unsigned int> group_start;	//< starting index for each group of nPhi sub-elements
};

#endif
