/* 
 * MixedSource.hh, part of the RotationShield program
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

#ifndef MIXEDSOURCE_HH
/// Make sure this header is only loaded once
#define MIXEDSOURCE_HH

#include "FieldSource.hh"
#include "LineSource.hh"
#include <vector>
using std::vector;

/// Combines several FieldSource objects into one (e.g. many LineSource wires into a Cos Theta coil)
class MixedSource: public FieldSource {
public:
    /// Constructor
    MixedSource(const string& nm = "MixedSource"): FieldSource(nm), displayColor(vec3(1.0,0.0,0.0)) {}
    
    /// Destructor
    virtual ~MixedSource() { clear(); }
    
    /// Field at a specified point
    vec3 fieldAt(const vec3& v) const;
    /// Field averaged over the given Line
    virtual vec3 fieldOverLine(Line l) const;
    /// Field averaged over the given Plane
    virtual vec3 fieldOverPlane(Plane p) const;
    /// Add a new FieldSource to the source currents
    virtual void addsource(const FieldSource* fs);
    /// Add FieldSources from another MixedSource
    void addsources(const MixedSource* MS);
    /// Remove all sources
    virtual void clear();
    /// Add FieldSource's specified in a file
    void loadSourcesFile(FILE* f, double scale);
    /// Get number of sources
    unsigned int nSources() const { return sources.size(); }
    
    /// Net current of the arrangement
    vec3 net_current() const { vec3 d = vec3(); for(unsigned int i=0; i<sources.size(); i++) d += sources[i]->net_current(); return d; }
    
    /// Add an arc of current (approximated by many LineSource segments)
    void arc(vec3 start, vec3 end, double j, int nsegs = 1);
    /// Add a loop of current (approximated by many LineSource segments)
    void loop(double z, double r, double j, int nsides = 32);
    
    /// Print info to stdout
    void display() const { printf("Multisource:\n"); for(unsigned int i=0; i<sources.size(); i++) sources[i]->display(); }
    /// Visualize the field source
    virtual void _visualize() const {
        vsr::setColor(displayColor[0],displayColor[1],displayColor[2]);
        for(unsigned int i=0; i<sources.size(); i++) sources[i]->_visualize();
    }
    
    vec3 displayColor;
    
private:

    vector<const FieldSource*> sources;    ///< component sources
};

#endif
