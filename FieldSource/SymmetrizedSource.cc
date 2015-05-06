#include "SymmetrizedSource.hh"
#include <cassert>

SymmetrizedSource::SymmetrizedSource(FieldSource* fs, bool parity): mySource(fs) {
    assert(mySource);
    mySource->retain();
    mySymmetry.rotation = mySource->getSymmetry().rotation;
    mySymmetry.parity = parity? 1 : -1;
}

SymmetrizedSource::~SymmetrizedSource() {
    mySource->release();
}
