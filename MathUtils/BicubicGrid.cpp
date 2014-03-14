#include "BicubicGrid.hh"
#include <cmath>
#include <cassert>
#include <stdio.h>


double eval_cubic(double y, const double* d) {
	return ( -0.5*d[0]*(1-y)*(1-y)*y
			+d[1]*(1-y)*(1-y*(1.5*y-1))
			-d[2]*y*(-0.5*(1-y)*(1-y)+y*(2*y-3))
			-0.5*d[3]*(1-y)*y*y );
}

double eval_cubic_deriv(double y, const double* d) {
	return ( -d[0]*(1.5*y*y -2*y +0.5)
			+d[1]*(4.5*y*y -5.*y)
			+d[3]*(1.5*y*y -y)
			-d[2]*(4.5*y*y -4*y -0.5));
}

//--------------------------------------------------------------------------

CubicGrid::CubicGrid(unsigned int nx): NX(nx) {
	setUserRange(0,1);
	data = new double[NX+4];
	bc = IB_LINEAR;
}

CubicGrid::~CubicGrid() {
	delete[] data;
}
	
double CubicGrid::operator()(double x) const {
	return _eval(sx*(x-ox));
}

double CubicGrid::deriv(double x) const {
	return _deriv(sx*(x-ox))*sx;
}

void CubicGrid::setUserRange(double r0, double r1, double e) {
	sx = (NX-1+2*e)/(r1-r0);
	ox = ((NX+1+e)*r0-(2-e)*r1)/(NX+2*e-1);
}

void CubicGrid::set(unsigned int x, double v) {
	assert(x<NX);
	_set(x+2,v);
}

double CubicGrid::_eval(double x) const {
	int ix = int(x);
	double fx = x-ix;
	
	// bounds check for out-of-range values
	if(ix<1 || ix > NX+1) {
		if(bc==IB_CYCLIC) ix = ((ix-2)+100*NX)%NX+2;
		else return 0;
	}
	
	assert(1<=ix && ix<=NX+1);
	return eval_cubic(fx, &data[ix-1]);
}

double CubicGrid::_deriv(double x) const {
	int ix = int(x);
	double fx = x-ix;
	
	// bounds check for out-of-range values
	if(ix<1 || ix > NX+1) {
		if(bc==IB_CYCLIC) ix = ((ix-2)+100*NX)%NX+2;
		else return 0;
	}
	
	assert(1<=ix && ix<=NX+1);
	return eval_cubic_deriv(fx, &data[ix-1]);
}

void CubicGrid::_set(unsigned int x, double v) {
	
	data[x] = v;
	
	// set boundary conditions
	
	if( 2<=x && x<=3 ) {
		if( bc==IB_LINEAR ) {
			data[1] = 2*data[2]-data[3];
			data[0] = 2*data[1]-data[2];
		}
		if( bc==IB_CYCLIC ) _set(NX+x,v);
	}
	if(x==2 && bc==IB_REPEAT) {  _set(1,v); _set(0,v); }
	
	if( NX <= x && x <= NX+1 ) {
		if( bc==IB_LINEAR ) {
			data[NX+2] = 2*data[NX+1]-data[NX];
			data[NX+3] = 2*data[NX+2]-data[NX+1];
		}
		if( bc==IB_CYCLIC ) _set(x-NX,v);
	}
	if(x==NX+1 && bc==IB_REPEAT) {  _set(NX+2,v); _set(NX+3,v); }
}


//--------------------------------------------------------------------------

BicubicGrid::BicubicGrid(unsigned int nx, unsigned int ny): NX(nx), NY(ny) {
	setUserRange(0,1,true);
	setUserRange(0,1,false);
	
	bc[0] = bc[1] = IB_REPEAT;
	
	data = new double*[NX+4];
	for(unsigned int i=0; i<NX+4; i++) {
		data[i] = new double[NY+4];
		for(unsigned int j=0; j<NY+4; j++) data[i][j]=0;
	}
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
	// sx*(r0-ox) = 2-e; sx*(r1-ox) = NX+1+e
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
	_set(x+2,y+2,v);
}

void BicubicGrid::_set(unsigned int x, unsigned int y, double v) {
	
	data[x][y] = v;
	
	// set boundary conditions
	
	if( 2<=x && x<=3 ) {
		if( bc[0]==IB_LINEAR ) {
			data[1][y] = 2*data[2][y]-data[3][y];
			data[0][y] = 2*data[1][y]-data[2][y];
		}
		if( bc[0]==IB_CYCLIC ) _set(NX+x,y,v);
	}
	if(x==2 && bc[0]==IB_REPEAT) {  _set(1,y,v); _set(0,y,v); }
	
	if( 2<=y &&  y<=3 ) {
		if( bc[1]==IB_LINEAR ) {
			data[x][1] = 2*data[x][2]-data[x][3];
			data[x][0] = 2*data[x][1]-data[x][2];
		}
		if( bc[1]==IB_CYCLIC ) _set(x, NY+y, v);
	}
	if(y==2 && bc[1]==IB_REPEAT) {  _set(x,1,v); _set(x,0,v); }
	
	if( NX <= x && x <= NX+1 ) {
		if( bc[0]==IB_LINEAR ) {
			data[NX+2][y] = 2*data[NX+1][y]-data[NX][y];
			data[NX+3][y] = 2*data[NX+2][y]-data[NX+1][y];
		}
		if( bc[0]==IB_CYCLIC ) _set(x-NX,y,v);
	}
	if(x==NX+1 && bc[0]==IB_REPEAT) {  _set(NX+2,y,v); _set(NX+3,y,v); }
	
	if( NY <= y && y <= NY+1 ) {
		if( bc[1]==IB_LINEAR ) {
			data[x][NY+2] = 2*data[x][NY+1]-data[x][NY];
			data[x][NY+3] = 2*data[x][NY+2]-data[x][NY+1];
		}
		if( bc[1]==IB_CYCLIC ) _set(x,y-NY,v);
	}
	if(y==NY+1 && bc[1]==IB_REPEAT) {  _set(x,NY+2,v); _set(x,NY+3,v); }
}

double BicubicGrid::eval_bicubic(double x, double y) const {
	double ypts[4];
	
	int ix = int(x);
	double fx = x-ix;
	int iy = int(y);
	double fy = y-iy;
	
	// bounds check for out-of-range values
	if(ix<1 || ix > NX+1) {
		if(bc[0]==IB_CYCLIC) ix = ((ix-2)+100*NX)%NX+2;
		else return 0;
	}
	if(iy<1 || iy > NY+1) {
		if(bc[1]==IB_CYCLIC) iy = ((iy-2)+100*NY)%NY+2;
		else return 0;
	}
	
	assert(1<=ix && ix<=NX+1 && 1<=iy && iy<=NY+1);
	
	for(int i=0; i<4; i++)
		ypts[i] = eval_cubic(fy, &data[ix-1+i][iy-1]);
	return eval_cubic(fx,ypts);
}

void BicubicGrid::printData() const {
	printf("----------- BicubicGrid %i x %i --------------\n",NX,NY);
	for(unsigned int nx=0; nx < NX+4; nx++) {
		for(unsigned int ny=0; ny<NY+4; ny++) {
			printf("%+0.3g\t",data[nx][ny]);
		}
		printf("\n");
	}
	printf("\n");
}

void BicubicGrid::minmax(double& mn, double& mx) const {
	mn = mx = data[2][2];
	for(unsigned int x=2; x<NX+2; x++) {
		for(unsigned int y=2; y<NY+2; y++) {
			double c = data[x][y];
			if(c<mn) mn = c;
			if(c>mx) mx = c;
		}
	}
	printf("Bicubic grid min=%g, max=%g\n",mn,mx);
}

void BicubicGrid::rescale(double a, double b) {
	for(unsigned int x=2; x<NX+2; x++)
		for(unsigned int y=2; y<NY+2; y++)
			_set(x,y,a+b*data[x][y]);
}

void BicubicGrid::scale_zrange(double a, double b) {
	double mn,mx;
	minmax(mn,mx);
	if(mx-mn > 1e-8)
		rescale( (a*mx-b*mn)/(mx-mn), (b-a)/(mx-mn) );
	else
		rescale(0.5*(a+b), 0);
}
