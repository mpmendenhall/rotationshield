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
	ReactiveSet(unsigned int nph=1): nPhi(nph), ixn_df(0) {}
	/// destructor
	virtual ~ReactiveSet() {}
	
	
	// Subclass me!
	//=====================================
	/// total number of degrees of freedom
	virtual unsigned int nDF() const = 0;
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0) = 0;
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip) { return false; }
	/// prepare incident state vector
	virtual void prepareIncident() { if(incidentState.size() != nDF()) incidentState = mvec(nDF()); }
	//=====================================
	
	/// clear states; start with only 1 active DF
	virtual void startInteractionScan() { setZeroState(); }
	/// set single degree of freedom to produce interactions from, clearing previous DF
	virtual void setInteractionDF(unsigned int DF, double v);
	/// set a DF without clearing previous DF
	virtual void setDF(unsigned int DF, double v) { finalState[DF] = v; _setDF(DF,v); }
	/// set zero state
	virtual void setZeroState() { setFinalState(mvec(nDF())); }
	/// set final state
	virtual void setFinalState(const mvec& v) { assert(v.size()==nDF()); finalState = v; _setDFv(v); }
	
	const unsigned int nPhi;	//< internal periodic symmetry
	mvec incidentState; 		//< non-interacting initial state vector
	mvec finalState;			//< final state after interactions
	
protected:
	
	// Subclass me!
	//=====================================
	/// called when a DF is set
	virtual void _setDF(unsigned int DF, double v) = 0;
	/// optional routine for setting entire state vector at once
	virtual void _setDFv(const mvec& v);
	//=====================================
	
	unsigned int ixn_df;				//< interaction DF currently set
};



//
//
//



/// linearly interacting DF grouped into logical sub-units; arranges shuffling indices into nPhi symmetry
class ReactiveUnitSet: public ReactiveSet {
public:
	/// constructor
	ReactiveUnitSet(unsigned int nph=1): ReactiveSet(nph) { group_start.push_back(0); }
	
	//=====================================
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return group_start.back(); }
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0);
	//=====================================
	
	
protected:
	
	// Subclass me!
	//=====================================
	/// set state for i^th sub-element
	virtual void setSubelDF(unsigned int el, unsigned int df, mdouble v) = 0;
	/// sub-element reaction: need this!
	virtual mvec subelReaction(unsigned int el, ReactiveSet* R) = 0;
	//=====================================
	
	
	/// additional routines for setting a DF value
	virtual void _setDF(unsigned int DF, double v=1.0);
	/// get DF index for sub-element DF
	unsigned int df_subindex(unsigned int el, unsigned int df) const { return group_start[el/nPhi] + df*nPhi + (el%nPhi); }
	/// get sub-element DF for global DF index
	void invert_index(unsigned int DF, unsigned int& el, unsigned int& df) const { el = df_subel[DF]; df = df_subel_df[DF]; }
	/// add a group of nPhi sub-elements with given number of DF
	void add_DF_group(unsigned int N);
	
	// interaction calculator cache variables
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



//
//
//



/// combine several ReactiveSets into one
class ReactiveSetCombiner: public ReactiveSet {
public:
	/// constructor
	ReactiveSetCombiner(unsigned int nph=1): ReactiveSet(nph) { set_cum_df.push_back(0); }
	/// append a new ReactiveSet
	virtual void addSet(ReactiveSet* R);
	
	
	//=====================================
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return set_cum_df.back(); }
	/// get DF for given phi reacting to state R
	virtual mvec getReactionTo(ReactiveSet* R, unsigned int phi = 0);
	/// respond to interaction protocol; return whether protocol recognized
	virtual bool queryInteraction(void* ip) { return mySets[ixn_set]->queryInteraction(ip); }
	/// prepare incident state vector
	virtual void prepareIncident();
	/// set single degree of freedom to produce interactions from, clearing previous DF
	virtual void setInteractionDF(unsigned int DF, double v);
	//=====================================
	
	/// reset interaction term counter
	virtual void startInteractionScan();
	
protected:

	//=====================================
	/// called when a DF is set
	virtual void _setDF(unsigned int DF, double v);
	/// optional routine for setting entire state vector at once
	virtual void _setDFv(const mvec& v);
	//=====================================
	
	unsigned int ixn_set;					//< set currently responsible for interacting
	std::vector<ReactiveSet*> mySets;		//< sub-units this class combined
	std::vector<unsigned int> set_cum_df;	//< cumulative DF of each subset
	std::vector<unsigned int> df_set;		//< which set each DF belongs to
};

#endif
