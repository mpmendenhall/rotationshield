/// \file ReactiveSet.hh \brief Virtual base class for collections of interacting degrees of freedom

#ifndef REACTIVESET_HH
/// Makes sure to only load this file once
#define REACTIVESET_HH 1

#include "Typedefs.hh"
#include <vector>
#include <cassert>

/// base class for a set of degrees of freedom
class ReactiveSet {
public:
	/// constructor
	ReactiveSet(unsigned int nph=1): nPhi(nph), ixn_df(0), ixn_ptcl(NULL) {}
	/// destructor
	virtual ~ReactiveSet() {}
	
	/// total number of degrees of freedom
	virtual unsigned int nDF() const = 0;
	
	/// set interaction protocol
	virtual bool set_protocol(void* ip) { ixn_ptcl = ip; return false; }
	/// reset interaction term counter
	virtual void startInteractionScan() { setZeroState(); setInteractionDF(0); }
	/// possibly arbitrarily ordered interaction matrix entries
	virtual mdouble nextInteractionTerm(unsigned int& i, unsigned int& j) = 0;
	/// set degree of freedom to produce interactions from
	virtual void setInteractionDF(unsigned int DF, double v=1.0);
	/// respond to interaction protocol
	virtual void queryInteraction() { }
	
	
	const unsigned int nPhi;	//< internal periodic symmetry

	mvec incidentState; 	//< non-interacting initial state vector
	
	/// set (final) state
	virtual void setFinalState(const mvec& v) { assert(v.size()==nDF()); finalState = v; }
	/// set zero state
	virtual void setZeroState() { setFinalState(mvec(nDF())); }

protected:

	mvec finalState;					//< final state after interactions
	mutable unsigned int rterm;			//< current interaction term
	unsigned int ixn_df;				//< interaction DF being queried
	void* ixn_ptcl;						//< current interaction protocol
	
};

/// linearly interacting DF grouped into logical sub-units; arranges shuffling indices into nPhi symmetry
class ReactiveUnitSet: public ReactiveSet {
public:
	/// constructor
	ReactiveUnitSet(unsigned int nph=1): ReactiveSet(nph) { group_start.push_back(0); }
	
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return group_start.back(); }
	
	/// reset interaction term counter
	virtual void startInteractionScan() { ReactiveSet::startInteractionScan(); ic_i = ic_di = 0; ic_v = subelReaction(); }
	/// possibly arbitrarily ordered interaction matrix entries
	virtual mdouble nextInteractionTerm(unsigned int& i, unsigned int& j);
	/// set degree of freedom to produce interactions from
	virtual void setInteractionDF(unsigned int DF, double v=1.0);
	
	/// set (final) state
	virtual void setFinalState(const mvec& v);
	
protected:
	
	/// get DF index for sub-element DF
	unsigned int df_subindex(unsigned int el, unsigned int df) const { return group_start[el/nPhi] + df*nPhi + (el%nPhi); }
	/// get sub-element DF for global DF index
	void invert_index(unsigned int DF, unsigned int& el, unsigned int& df) const { el = df_subel[DF]; df = df_subel_df[DF]; }
	
	/// add a group of nPhi sub-elements with given number of DF
	void add_DF_group(unsigned int N);
	
	/// set state for i^th sub-element
	virtual void setSubelDF(unsigned int el, unsigned int df, mdouble v) = 0;
		
	/// sub-element reaction to RS via protocol: need this!
	virtual mvec subelReaction() = 0;
	
	// interaction calculator cache variables
	mvec ic_v;				//< cached interaction matrix between sub-elements
	unsigned int ic_i;		//< interaction scan element being examined
	unsigned int ic_di;		//< interaction scan element DF being examined
	unsigned int ixn_el;	//< sub-element currently set for interaction study
	unsigned int ixn_el_df;	//< sub-element's DF currently set for interaction study
	
	// groups
	unsigned int n_subels() const { return group_DF.size()*nPhi; }
	std::vector<unsigned int> group_DF;		//< number of DF for elements in this group
	std::vector<unsigned int> group_start;	//< starting index for each group of nPhi sub-elements
	
	// individual DF
	std::vector<unsigned int> df_subel;		//< sub-element number corresponding to each DF
	std::vector<unsigned int> df_subel_df;	//< sub-element's DF corresponding to each DF
};

// TODO: combining RS

#endif
