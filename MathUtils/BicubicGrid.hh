#ifndef BICUBICGRID_HH
#define BICUBICGRID_HH 1

/// simple bicubic grid; lacking flexibility, but faster than general interpolator schemes
class BicubicGrid {
public:
	/// constructor
	BicubicGrid(unsigned int nx, unsigned int ny);
	/// destructor
	~BicubicGrid();
	
	const unsigned int NX;
	const unsigned int NY;
	
	/// evaluate at (x,y) in user coordinates
	double operator()(double x, double y) const;
	/// set scale factors for user range, points at ends
	void setUserRange(double r0, double r1, bool xdirection, double e = 0);
	/// boundary conditions for interpolation
	enum IBC {
		IB_CYCLIC,	//< cyclic edges
		IB_ZERO,	//< zero-pad edges
		IB_LINEAR	//< linear approach to edges
	} bc[2];	//< for each axis
	/// set value at point x,y
	void set(unsigned int x, unsigned int y, double v);
	
protected:

	/// eval cubic for y in 0,1, points at -1,0,1,2
	double eval_cubic(double y, double* d) const;
	/// eval bicubic, x & y scaled to data range
	double eval_bicubic(double x, double y) const;
	
	double sx,sy,ox,oy;		//< user coordinate locations of first and last point in each dimension
	double** data;			//< data with edge guard values
};

#endif
