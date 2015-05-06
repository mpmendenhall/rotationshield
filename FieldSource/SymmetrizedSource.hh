#ifndef SYMMETRIZEDSOURCE_HH
#define SYMMETRIZEDSOURCE_HH

#include "FieldSource.hh"

/// Symmetrized/antisymmetrized source configuration f(x) +- f(-x)
class SymmetrizedSource: public FieldSource {
public:
    /// Constructor
    SymmetrizedSource(FieldSource* fs, bool parity);
    /// Destructor
    virtual ~SymmetrizedSource();
    
    /// (symmetrized) field at given point
    virtual vec3 fieldAt(const vec3& v) const { return mySource->fieldAt(v) + mySource->fieldAt(-v)*mySymmetry.parity; }
    /// Calculate effective net current of field source
    virtual vec3 net_current() const { return mySymmetry.parity > 0? mySource->net_current()*2 : vec3(); }
    /// Print info to stdout
    virtual void display() const { printf("[Symmetrized %+i]",mySymmetry.parity); mySource->display(); }
    /// Visualize source (TODO with symmetrized copy)
    virtual void _visualize() const { mySource->_visualize(); }
    
protected:
    FieldSource* mySource;
};

#endif
