/* 
 * InterpolatorHelper.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "InterpolationHelper.hh"
#include <algorithm>
#include <cassert>

InterpolationHelper::InterpolationHelper(): DataSequence(BC_INFINITE), myInterpolator(NULL) {
    setInterpolatorMethod(&CubiTerpolator::newCubiTerpolator);
    cum_sum_dpts.push_back(0);
}

InterpolationHelper::~InterpolationHelper() {
    clearSubinterps();
    if(myInterpolator) delete myInterpolator;
}

double InterpolationHelper::eval(const double* x) const {
    assert(x);
    return myInterpolator->eval(x);
}

void InterpolationHelper::setInterpolatorMethod(Interpolator* (*makeInterp)(DataSequence*), unsigned int nDeep) {
    if(nDeep) {
        for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++)
            (*it)->setInterpolatorMethod(makeInterp,nDeep-1);
    } else {
        Interpolator* newInterp = makeInterp(this);
        if(myInterpolator) {
            newInterp->scale = myInterpolator->scale;
            newInterp->offset = myInterpolator->offset;
            delete myInterpolator;
        }
        myInterpolator = newInterp;
    }
}

void InterpolationHelper::setBoundaryCondition(BoundaryCondition b, unsigned int nDeep) {
    if(!nDeep) bc = b;
    else {
        for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++)
            (*it)->setBoundaryCondition(b, nDeep-1);
    }
}

void InterpolationHelper::setupDataGrid(const unsigned int* n0, const unsigned int* n1) {
    assert(n1>n0);
    clearSubinterps();
    myData.clear();
    if(n1 == n0+1) {
        myData = std::vector<double>(*n0);
    } else {
        for(unsigned int i = 0; i < *n0; i++) {
            subInterpolators.push_back(new InterpolationHelper());
            subInterpolators.back()->setupDataGrid(n0+1,n1);
        }
    }
    recalc_structure(false);
}

double& InterpolationHelper::operator[](const unsigned int* n) {
    assert(n);
    if(isBottom()) {
        assert(*n < myData.size());
        return myData[*n];
    } else {
        assert(*n < subInterpolators.size());
        return (*subInterpolators[*n])[n+1];
    }
}

double& InterpolationHelper::operator[](unsigned int n) {
    if(isBottom()) {
        assert(n < myData.size());
        return myData[n];
    } else {
        unsigned int m = from_flat_index(n);
        return (*subInterpolators[m])[n];
    }
}

void InterpolationHelper::zeroData() {
    if(isBottom()) std::fill(myData.begin(), myData.end(), 0);
    else {
        for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++)
            (*it)->zeroData();
    }
}

void InterpolationHelper::setData(const double* x0) {
    if(isBottom())
         std::copy(x0, x0+myData.size(), myData.begin());
    else {
        for(unsigned int i=0; i<subInterpolators.size(); i++)
            subInterpolators[i]->setData(x0+cum_sum_dpts[i]);
    }
}

const InterpolationHelper& InterpolationHelper::getSubHelper(unsigned int nDeep, const unsigned int* n) const {
    if(!nDeep) return *this;
    assert(n);
    assert(n[0] < subInterpolators.size());
    return subInterpolators[n[0]]->getSubHelper(nDeep-1,n+1);
}

std::vector<InterpolationHelper*> InterpolationHelper::getSubHelpers(unsigned int nDeep) {
    std::vector<InterpolationHelper*> v;
    if(!nDeep) v.push_back(this);
    else {
        for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++) {
            std::vector<InterpolationHelper*> v2 = (*it)->getSubHelpers(nDeep-1);
            v.insert(v.end(), v2.begin(), v2.end());
        }
    }
    return v;
}


double InterpolationHelper::valueAt(int i, void* xopts) const {
    unsigned int ci = coerce(i);
    if(ci < myData.size()) return myData[ci];
    if(ci < subInterpolators.size()) return subInterpolators[ci]->eval((double*)xopts);
    assert(false);
    return 0;
}

void InterpolationHelper::clearSubinterps() {
    for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++)
        delete *it;
    subInterpolators.clear();
}

void InterpolationHelper::recalc_structure(bool recurse) {
    cum_sum_dpts.clear();
    cum_sum_dpts.push_back(0);
    if(isBottom()) {
        npts = myData.size();
        cum_sum_dpts.push_back(myData.size());
    } else {
        npts = subInterpolators.size();
        for(std::vector<InterpolationHelper*>::iterator it = subInterpolators.begin(); it != subInterpolators.end(); it++) {
            if(recurse) (*it)->recalc_structure(true);
            cum_sum_dpts.push_back(cum_sum_dpts.back()+(*it)->n_pts());
        }
    }
}

unsigned int InterpolationHelper::from_flat_index(unsigned int& n) const {
    std::vector<unsigned int>::const_iterator it = std::upper_bound(cum_sum_dpts.begin(), cum_sum_dpts.end(), n);
    assert(it != cum_sum_dpts.end());
    unsigned int m = it-cum_sum_dpts.begin()-1;
    assert(m<cum_sum_dpts.size());
    n -= cum_sum_dpts[m];
    return m;
}



