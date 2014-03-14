#ifndef BICUBICGRID_HH
#define BICUBICGRID_HH 1

/// evaluate cubic spline for y in 0,1, points at -1,0,1,2
double eval_cubic(double y, const double* d);
/// derivative of cubic spline
double eval_cubic_deriv(double y, const double* d);

/// boundary conditions for interpolation
enum IBC {
		IB_CYCLIC,	//< cyclic edges
		IB_ZERO,	//< zero-pad edges
		IB_LINEAR,	//< linear approach to edges
		IB_REPEAT	//< repeat end value
};

/// simple 1D bicubic interpolator; faster than general interpolator schemes
class CubicGrid {
public:
	/// constructor
	CubicGrid(unsigned int nx);
	/// destructor
	~CubicGrid();
	
	const unsigned int NX;	//< number of grid points
	IBC bc;					//< boundary condition
	
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
	
	double sx,ox;	//< user coordinate locations of first and last point in each dimension
	double* data;	//< data with edge guard values
};


/// simple bicubic grid; faster than general interpolator schemes
class BicubicGrid {
public:
	/// constructor
	BicubicGrid(unsigned int nx, unsigned int ny);
	/// destructor
	~BicubicGrid();
	
	const unsigned int NX;	//< number of grid points in x direction
	const unsigned int NY;	//< number of grid points in y direction
	
	/// evaluate at (x,y) in user coordinates
	double operator()(double x, double y) const;
	
	/// set scale factors for user range, points at ends
	void setUserRange(double r0, double r1, bool xdirection, double e = 0);
	
	IBC bc[2];	//< boundary condition for each axis
	
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
	/// set value at point x,y, in internal coordinates
	void _set(unsigned int x, unsigned int y, double v);
	
	double sx,sy,ox,oy;		//< user coordinate locations of first and last point in each dimension
	double** data;			//< data with edge guard values
};

#endif
