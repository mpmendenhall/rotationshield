#ifndef SYMMETRYSPEC_HH
#define SYMMETRYSPEC_HH

/// Specification for symmetries that a field/surface can have
class SymmetrySpec {
public:
    /// Constructor
    SymmetrySpec() { }
    
    /// combined symmetry
    SymmetrySpec& operator+=(const SymmetrySpec& rhs) { rotation &= rhs.rotation; parity &= (parity == rhs.parity); return *this; } 
    
    bool rotation = false;      ///< symmetry under rotation around z
    int parity = 0;             ///< symmetry under (x,y,z)->(-x,-y,-z); -1, 0, or 1
};

#endif
