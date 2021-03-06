/* 
 * BicubicGrid.hh, part of the RotationShield program
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

#ifndef BICUBICGRID_HH
/// Make sure this header is only loaded once
#define BICUBICGRID_HH

/// evaluate cubic spline for y in 0,1, points at -1,0,1,2
double eval_cubic(double y, const double* d);
/// derivative of cubic spline
double eval_cubic_deriv(double y, const double* d);

/// boundary conditions for interpolation
enum IBC {
        IB_CYCLIC,      ///< cyclic edges
        IB_ZERO,        ///< zero-pad edges
        IB_LINEAR,      ///< linear approach to edges
        IB_REPEAT       ///< repeat end value
};

/// simple 1D bicubic interpolator; faster than general interpolator schemes
class CubicGrid {
public:
    /// constructor
    CubicGrid(unsigned int nx);
    /// destructor
    ~CubicGrid();
    
    const unsigned int NX;      ///< number of grid points
    IBC bc;                     ///< boundary condition
    
    /// evaluate at x in user coordinates
    double operator()(double x) const;
    /// evaluate derivative at x in user coordinates
    double deriv(double x) const;
    
    /// set scale factors for user range, points at ends
    void setUserRange(double r0, double r1, double e = 0);
    /// set value at point x
    void set(unsigned int x, double v);
    
protected:

    /// eval cubic, x scaled to data range
    double _eval(double x) const;
    /// eval derivative, x scaled to data range
    double _deriv(double x) const;
    
    /// set value at point x,y, in internal coordinates
    void _set(unsigned int x, double v);
    
    double sx,ox;       ///< user coordinate locations of first and last point in each dimension
    double* data;       ///< data with edge guard values
};


/// simple bicubic grid; faster than general interpolator schemes
class BicubicGrid {
public:
    /// constructor
    BicubicGrid(unsigned int nx, unsigned int ny);
    /// destructor
    ~BicubicGrid();
    
    const unsigned int NX;      ///< number of grid points in x direction
    const unsigned int NY;      ///< number of grid points in y direction
    
    /// evaluate at (x,y) in user coordinates
    double operator()(double x, double y) const;
    /// evaluate derivative along one axis in user coordinates
    double deriv(double x, double y, bool xdirection) const;
    
    /// set scale factors for user range, points at ends
    void setUserRange(double r0, double r1, bool xdirection, double e = 0);
    
    IBC bc[2];    ///< boundary condition for each axis
    
    /// set value at point x,y
    void set(unsigned int x, unsigned int y, double v);
    
    /// print data grid to screen
    void printData() const;
    
    /// get maximum and minimum value
    void minmax(double& mn, double& mx) const;
    /// apply a+bx linear data transform
    void rescale(double a, double b);
    /// scale to fit range
    void scale_zrange(double a, double b);
    
protected:

    /// eval bicubic, x & y scaled to data range
    double eval_bicubic(double x, double y) const;
    /// evaluate derivative, x & y scaled to data range
    double eval_deriv(double x, double y, bool xdirection) const;
    /// set value at point x,y, in internal coordinates
    void _set(unsigned int x, unsigned int y, double v);
    
    double sx,sy,ox,oy; ///< user coordinate locations of first and last point in each dimension
    double** data;      ///< data with edge guard values
};

#endif
