/* 
 * ScanRange.hh, part of the RotationShield program
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

/// \file "ScanRange.hh" \brief Class for scanning a variable over an interval
#ifndef SCANRANGE_HH
/// Make sure this header is only loaded once
#define SCANRANGE_HH

/// Class for scanning over an interval in a specified number of uniform steps
class ScanRange {
public:
    /// Constructor
    /**    \param start first point in scan
     \param end last point in scan
     \param npts number of points in scan */
    ScanRange(double start, double end, int npts): nmax(npts), s(start), e(end), n(0) { }
    
    /// Destructor (nothing to do)
    ~ScanRange() {}
    
    /// Step to next point in interval
    double next() { if(nmax>1) return s + (e-s)*double(n++)/double(nmax - 1); else return 0.5*(s+e)*(++n); }
    /// Get current location
    double current() const { if(nmax>1) return s + (e-s)*double(n-1)/double(nmax - 1); else return 0.5*(s+e); }
    /// Check whether there are more points to go after this one
    bool goOn() const { return n<=nmax; }
    /// Return the current step number #n
    int getN() const { return n; }
    /// Print the status of the scan to stdout
    void printStatus() const { printf("%g [%i/%i]",current(),n,nmax); }
    
    const int nmax; ///< total number of steps over interval
    const double s; ///< Starting point of interval
    const double e; ///< Ending point of interval
private:
    int n; ///< Step number in [0,#nmax)
};


#endif
