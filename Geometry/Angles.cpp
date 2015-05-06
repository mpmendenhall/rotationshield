/* 
 * Angles.cpp, part of the RotationShield program
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

#include "Angles.hh"
#include <cassert>

float normalizeAngle(float a, float theta0) {
    if(a<theta0) a += 2*M_PI*int(1+(theta0-a)/(2*M_PI));
    if(a>=theta0+2*M_PI) a -= 2*M_PI*int((a-theta0)/(2*M_PI));
    return a;
}

//-----------------------------------------

angular_interval::angular_interval(double t0, double t1, bool nrm): th0(t0), th1(t1) {
    if(nrm) normalize();
}

void angular_interval::normalize() {
    int n2pi = th0/(2*M_PI);
    if(th0<0) n2pi--;
    th0 -= n2pi*2*M_PI;
    th1 -= n2pi*2*M_PI;
}

void angular_interval::add(double a) {
    th0 += a;
    th1 += a;
    normalize();
}

void Angular_Interval_Set::subtract_interval(angular_interval i) {
    // do nothing for trivial intervals
    if(!endpts.size() || i.th0==i.th1) return;
    
    // special case: split wrap-around intervals
    if(i.th1 > 2*M_PI) {
        subtract_interval(angular_interval(i.th0, 2*M_PI));
        subtract_interval(angular_interval(0, i.th1-2*M_PI));
        return;
    }
        
    angular_interval ilo(i.th0,i.th0);
    angular_interval iup(i.th1,4*M_PI,false);
    std::set<angular_interval>::iterator it0 = endpts.lower_bound(ilo);
    std::set<angular_interval>::iterator it1 = endpts.upper_bound(iup);
    
    // identify partial intervals to delete and restore truncated
    std::set<angular_interval>::iterator it0m = it0;
    std::set<angular_interval>::iterator it1m = it1;
    bool tr0 = (it0 != endpts.begin() && (--it0m)->th1 > i.th0);
    bool tr1 = (it1 != endpts.begin() && (--it1m)->th1 > i.th1);
        
    std::vector<angular_interval> partials;
    if(tr0) {
        partials.push_back(angular_interval(it0m->th0,i.th0));
        --it0;
    }
    if(tr1) {
        partials.push_back(angular_interval(i.th1,it1m->th1));
    }
    
    endpts.erase(it0,it1);
    
    for(std::vector<angular_interval>::const_iterator it = partials.begin(); it != partials.end(); it++)
        endpts.insert(*it);
}

void Angular_Interval_Set::Angular_Interval_Set::add_interval(angular_interval i) {
    // do nothing for trivial intervals
    if(i.th0==i.th1) return;
    
    // special case: split wrap-around intervals
    if(i.th1 > 2*M_PI) {
        add_interval(angular_interval(i.th0,2*M_PI));
        add_interval(angular_interval(0,i.th1-2*M_PI));
        return;
    }
    
    angular_interval iup(i.th1,4*M_PI);
    std::set<angular_interval>::iterator it0 = endpts.lower_bound(i);
    std::set<angular_interval>::iterator it1 = endpts.upper_bound(iup);
    
    // intersection with existing intervals
    std::set<angular_interval>::iterator it0m = it0;
    std::set<angular_interval>::iterator it1m = it1;
    bool tr0 = (it0 != endpts.begin() && (--it0m)->th1 >= i.th0);
    bool tr1 = (it1 != endpts.begin() && (--it1m)->th1 >= i.th1);
    
    if(tr0) { i.th0 = it0m->th0; it0--; }
    if(tr1) i.th1 = it1m->th1;
    
    endpts.erase(it0,it1);
    
    //std::cout << "Adding " << i << std::endl;
    endpts.insert(i);
}

std::vector<angular_interval> Angular_Interval_Set::get_intervals(bool merge_wrap) const {
    std::vector<angular_interval> v;
    for(std::set<angular_interval>::const_iterator it = endpts.begin(); it != endpts.end(); it++)
        v.push_back(*it);
    if(merge_wrap && v.size()>=2 && v[0].th0 == 0 && v.back().th1 == 2*M_PI) {
        v[0].th0 = v.back().th0;
        v[0].th1 += 2*M_PI;
        v.pop_back();
    }
    return v;
}


//-----------------------------

std::ostream& operator<<(std::ostream& o, const angular_interval& i) {
    o << "[" << i.th0/M_PI << " pi, " << i.th1/M_PI << " pi]";
    return o;
}

std::ostream& operator<<(std::ostream& o, const Angular_Interval_Set& S) {
    std::vector<angular_interval> v = S.get_intervals();
    if(!v.size()) o << "{}";
    for(std::vector<angular_interval>::iterator it = v.begin(); it != v.end(); it++) {
        if(it != v.begin()) o << " U ";
        o << *it;
    }
    return o;
}
