#include "BicubicGrid.hh"
#include <cassert>

BicubicGrid::BicubicGrid(unsigned int nx, unsigned int ny): NX(nx), NY(ny) {
	setUserRange(0,1,true);
	setUserRange(0,1,false);
	
	bc[0] = bc[1] = IB_CYCLIC;
	
	data = new double*[NX];
	for(unsigned int i=0; i<NX; i++)
		data[i] = new double[NY];
}

BicubicGrid::~BicubicGrid() {
	for(unsigned int i=0; i<NX; i++)
		delete[] data[i];
	delete[] data;
}

double BicubicGrid::operator()(double x, double y) const {
	return eval_bicubic(sx*(x-ox), sy*(y-oy));
}
	
void BicubicGrid::setUserRange(double r0, double r1, bool xdirection, double e) {
	// sx*(r0-ox) = (2-e); sx*(r1-ox) = NX+1+e
	if(xdirection) {
		sx = (NX-1+2*e)/(r1-r0);
		ox = ((NX+1+e)*r0-(2-e)*r1)/(NX+2*e-1);
	} else {
		sy = (NY-1+2*e)/(r1-r0);
		oy = ((NY+1+e)*r0-(2-e)*r1)/(NY+2*e-1);
	}
}

void BicubicGrid::set(unsigned int x, unsigned int y, double v) {
	assert(x<NX && y<NY);
	
	x += 2;
	y += 2;
	data[x][y] = v;
	
	// set boundary conditions
	
	if( x<=3 && bc[0]==IB_LINEAR ) {
		data[1][y] = 2*data[2][y]-data[3][y];
		data[0][y] = 2*data[1][y]-data[2][y];
	}
	if( x>=NX && bc[0]==IB_LINEAR ) {
		data[NX+2][y] = 2*data[NX+1][y]-data[NX][y];
		data[NX+3][y] = 2*data[NX+2][y]-data[NX+1][y];
	}
	if( y<=3 && bc[1]==IB_LINEAR ) {
		data[x][1] = 2*data[x][2]-data[x][3];
		data[x][0] = 2*data[x][1]-data[x][2];
	}
	if( y>=NY && bc[1]==IB_LINEAR ) {
		data[x][NY+2] = 2*data[x][NY+1]-data[x][NY];
		data[x][NY+3] = 2*data[x][NY+2]-data[x][NY+1];
	}
	
	if( x<=3 && bc[0]==IB_CYCLIC ) data[NX+x][y] = data[x][y];
	if( y<=3 && bc[1]==IB_CYCLIC ) data[x][NY+y] = data[x][y];
	if( x>=NX && bc[0]==IB_CYCLIC ) data[x-NX][y] = data[x][y];
	if( y>=NY && bc[1]==IB_CYCLIC ) data[x][y-NY] = data[x][y];
	
}

double BicubicGrid::eval_cubic(double y, double* d) const {
	return ( -0.5*d[0]*(1-y)*(1-y)*y
			+d[1]*(1-y)*(1-y*(1.5*y-1))
			-d[2]*y*(-0.5*(1-y)*(1-y)+y*(2*y-3))
			-0.5*d[2]*(1-y)*y*y );
}

double BicubicGrid::eval_bicubic(double x, double y) const {
	double ypts[4];
	
	int ix = int(x);
	double fx = x-ix;
	int iy = int(y);
	double fy = y-iy;
	
	if(ix<1 || iy<1 || ix > NX+1 || iy > NY+1) return 0; // bounds check for out-of-range values
	
	for(int i=0; i<4; i++)
		ypts[i] = eval_cubic(fy, &data[ix-1+i][iy-1]);
	return eval_cubic(fx,ypts);
}
