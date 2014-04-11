/// \file "interpolator.hh" \brief Adaptive mesh interpolator

#ifndef INTERPOLATOR_HH
/// Make sure this header is only loaded once
#define INTERPOLATOR_HH

#include "Visr.hh"
#include <math.h>

template<class V>
class TopMesh;

template<class V>
class InterpoBasis
{
public:
	InterpoBasis(int ns, double (*sf)(double)): nsupport(ns*ns), splinefunc(sf),
	supportx((int*)malloc(ns*ns*sizeof(int))), supporty((int*)malloc(ns*ns*sizeof(int)))
	{
		for(int i=0; i<ns; i++)
		{
			for(int j=0; j<ns; j++)
			{
				supportx[ns*i+j] = 1+i-ns/2;
				supporty[ns*i+j] = 1+j-ns/2;
			}
		}		
	}
	virtual ~InterpoBasis() { free(supportx); free(supporty); }
	
	virtual V eval(double x, double y, V** coeffs) const
	{
		V a;
		for(int i=0; i<nsupport; i++)
		{
			a = a + (*coeffs[i])*splinefunc(double(supportx[i])-x-0.5)*splinefunc(double(supporty[i])-y-0.5);
		}
		return a;
	}
	double (*splinefunc)(double);
	int nsupport;
	int* supportx;
	int* supporty;
};


double linInterp(double x)
{
	x = fabs(x);
	if(x>1) return 0;
	return 1-x;
}
	
double cubInterp(double x)
{
	if(x<0) x=-x;
	if(x>2) return 0;
	if(x>1) return ( ( -1.0/3.0 * (x-1) + 4.0/5.0     ) * (x-1) -   7.0/15.0 ) * (x-1) ;
	return ( ( x - 9.0/5.0 ) * x -   1.0/5.0     ) * x + 1.0;
}

double sp16(double x)
{
	if(x<0) x=-x;
	if(x>2) return 0;
	double A = -0.75;
	if(x>1) return (( A * x - 5.0 * A ) * x + 8.0 * A ) * x - 4.0 * A;
	return (( A + 2.0 )*x - ( A + 3.0 ))*x*x +1.0;	
}

double sp36(double x)
{
	if(x<0) x=-x;
	if(x>3) return 0;
	if(x>2) return ( (    1.0/11.0  * (x-2) -  45.0/ 209.0 ) * (x-2) +  26.0/ 209.0  ) *(x-2);
	if(x>1) return ( ( -  6.0/11.0  * (x-1) + 270.0/ 209.0 ) * (x-1) - 156.0/ 209.0  ) *(x-1);
	return ( (   13.0/11.0  * x - 453.0/ 209.0 ) * x -   3.0/ 209.0  ) * x + 1.0;
}

InterpoBasis<vec3>* spline36basis = new InterpoBasis<vec3>(6,&sp36);








template<class V>
class dataTree
{
public:
	dataTree(): subTrees(NULL), val(NULL) {}
	dataTree(V* v): val(v), subTrees(NULL) {}
	dataTree(V v): val((V*)malloc(sizeof(V))), subTrees(NULL) { *val = v; }
	
	virtual ~dataTree()
	{
		if(subTrees)
		{
			for(int i=0; i<4; i++) delete(subTrees[i]);
			free(subTrees);
		}
		
		if(val) free(val);
	}
	
	void branch()
	{
		subTrees = (dataTree<V>**)malloc(4*sizeof(V*));
		for(int i=1; i<4; i++) subTrees[i] = new dataTree();
		subTrees[0] = new dataTree(val);
	}
	
	virtual V* dtGetVal(const int d, const int lx, const int ly)
	{
		if(!d) return val;
		if(!subTrees) return NULL;
		int c = (lx/(1<<(d-1)))+2*(ly/(1<<(d-1)));
		return subTrees[c]->dtGetVal(d-1,lx%(1<<(d-1)),ly%(1<<(d-1)));
	}
	
	virtual V* setVal(const int d, const int lx, const int ly, V v)
	{
		if(!d) { val = (V*)malloc(sizeof(V)); (*val)=v; return val; }
		if(!subTrees) branch();
		int c = (lx/(1<<(d-1)))+2*(ly/(1<<(d-1)));
		return subTrees[c]->setVal(d-1,lx%(1<<(d-1)),ly%(1<<(d-1)),v);
	}
	
	virtual void forgetVal() { val = NULL; }
	
	// placeholders for InterpSquare methods
	virtual void checkSubdivide(int d) {}
	virtual int doSplit(const int d) { return 0; }
	virtual V interpolate(const int d, const int lx, const int ly) { return V(); }
	virtual void visualize(bool istop = true) {}
	virtual void visualize_cyl(double r = 0) {}
	virtual void markForSplit(int d, int lx, int ly) {}
	
protected:
	dataTree<V>** subTrees;
	V* val;
};





template<class V>
class InterpSquare: public dataTree<V>
{
protected:
	TopMesh<V>* tm;
	V** coeffs;
	bool split[4];
	const int mydepth, absx, absy;
	
public:
	/// Constructor
	InterpSquare(TopMesh<V>* TM, const int d, const int ax, const int ay, V* v):
	dataTree<V>(v), tm(TM), coeffs((V**)malloc(TM->b->nsupport*sizeof(V*))),
	mydepth(d), absx(ax), absy(ay)
	{
		for(int i=0; i<4; i++)
			split[i] = false;
	}
	
	/// Destructor
	virtual ~InterpSquare() { }
	
	virtual V interpolate(double lx, double ly) const
	{ 
		int c = (lx<0)+2*(ly<0);
		if(split[c])
		{
			double nx = 2*lx-0.5; if(nx>=0.5) nx -= 1.0;
			double ny = 2*ly-0.5; if(nx>=0.5) ny -= 1.0;
			return ((InterpSquare*)this->subTrees[c])->interpolate(nx,ny);
		}
		return tm->b->eval(lx,ly,coeffs);
	}
	
	virtual V interpolate(const int d, const int lx, const int ly)
	{
		if(!d) return *(this->dtGetVal(0,0,0));
		int c = (lx/(1<<(d-1)))+2*(ly/(1<<(d-1)));
		if(split[c]) return this->subTrees[c]->interpolate(d-1,lx%(1<<(d-1)),ly%(1<<(d-1)));
		return shallowInterpolate(d,lx,ly);
	}
		
	virtual inline V shallowInterpolate(const int d, const int lx, const int ly) const
	{
		return tm->b->eval(location(d,lx),location(d,ly),coeffs);
	}
	
	inline double location(const int d, const int lx) const
	{
		return double(lx)/double(1<<d)-0.5;
	}
	
	inline double absLocation(const int d, const int n, const bool xaxis) const
	{
		return tm->absLocation(mydepth+d, absx*(1<<d)+n, xaxis);
	}
	
	virtual V* getVal(const int d, const int lx, const int ly)
	{
		V* v = ((dataTree<V>*)this)->dtGetVal(d,lx,ly);
		if(v) return v;
		return setVal(d,lx,ly,interpolate(d,lx,ly));
	}
		
	virtual void checkSubdivide(int d)
	{
		if(d)
		{
			for(int i=0; i<4; i++)
				if(split[i])
					this->subTrees[i]->checkSubdivide(d-1);
			return;
		}
		
		for(int dx=0; dx<2; dx++)
		{
			for(int dy=0; dy<2; dy++)
			{
				if(!dx && !dy) continue;
				V actualval = tm->evalAt(mydepth+1,2*absx+dx,2*absy+dy);
				setVal(1,dx,dy,actualval);
				V interpval = shallowInterpolate(1,dx,dy);
				if(!(tm->closeEnough(actualval,interpval)))
					tm->markSplitAffected(mydepth+1,2*absx+dx,2*absy+dy);
			}
		}
	}
	
	virtual int doSplit(const int d)
	{
		int nsplit=0;
		
		if(!d)
		{
			for(int i = 0; i<4; i++)
			{
				if(split[i])
				{
					InterpSquare* s = new InterpSquare(tm, mydepth+1, 2*absx+i%2, 2*absy+i/2,this->dtGetVal(1,i%2,i/2));
					this->subTrees[i]->forgetVal();
					delete(this->subTrees[i]);
					this->subTrees[i] = s;
					nsplit++;
				}
			}
			return nsplit;
		}
		
		for(int i = 0; i<4; i++)
			nsplit += this->subTrees[i]->doSplit(d-1);
		return nsplit;
	}
		
	virtual void markForSplit(int d, int lx, int ly)
	{
		if(d==1)
			split[lx+2*ly] = true;
		else
		{
			int c = (lx/(1<<(d-1)))+2*(ly/(1<<(d-1)));
			if(split[c]) this->subTrees[c]->markForSplit(d-1,lx%(1<<(d-1)),ly%(1<<(d-1)));
		}
	}
	
	virtual void display(int v = 0)
	{
		printf("InterpSquare %i %i,%i\n",mydepth,absx,absy);
		if(v>5)
		{
			for(int i=0; i<tm->b->nsupport; i++)
				coeffs[i]->display();
		}
	}
	
	virtual void _visualize() const
	{
		
		int nsub = 0;
		for(int i=0; i<4; i++)
		{
			if(split[i])
			{
				this->subTrees[i]->_visualize();
				nsub++;
			}
		}
		
		if(!nsub)
		{
			vsr::setColor(0.0,0.0,1.0);
			float xyzcorners[12];
			for(int dx=0; dx<2; dx++)
			{
				for(int dy=0; dy<2; dy++)
				{
					xyzcorners[3*(dx+2*dy)+0] = tm->absLocation(mydepth,absx+dx,true);
					xyzcorners[3*(dx+2*dy)+1] = tm->absLocation(mydepth,absy+dy,false);
					xyzcorners[3*(dx+2*dy)+2] = tm->getVal(mydepth,absx+dx,absy+dy)->mag();
				}
			}
			vsr::quad(xyzcorners);
		
		}
	}
	
	virtual void visualize_cyl(double r = 0.75)
	{
		
		int nsub = 0;
		for(int i=0; i<4; i++)
		{
			if(split[i])
			{
				this->subTrees[i]->visualize_cyl(r);
				nsub++;
			}
		}
		
		if(!nsub)
		{
			
			float xyzcorners[12];
			double x,y,z;
			for(int i=0; i<3; i++)
			{
				if(i==0)
					vsr::setColor(1.0,0.0,0.0);
				else if(i==1)
					vsr::setColor(0.0,1.0,0.0);
				else
					vsr::setColor(0.0,0.0,1.0);
				
				for(int dx=0; dx<2; dx++)
				{
					for(int dy=0; dy<2; dy++)
					{
						x = tm->absLocation(mydepth,absx+dx,true);
						y = tm->absLocation(mydepth,absy+dy,false);
						z = r+(*(tm->getVal(mydepth,absx+dx,absy+dy)))[i];
						xyzcorners[3*(dx+2*dy)+0] = z*cos(x);
						xyzcorners[3*(dx+2*dy)+1] = z*sin(x);
						xyzcorners[3*(dx+2*dy)+2] = y;
					}
				}
				vsr::quad(xyzcorners);
			}
		}
	
	}
	
	virtual void updateCoeffs(int d)
	{
		if(!d)
		{
			for(int i=0; i<tm->b->nsupport; i++)
				coeffs[i] = tm->getVal(mydepth,absx+tm->b->supportx[i],absy+tm->b->supporty[i]);
		}
		else
		{
			for(int i=0; i<4; i++)
				if(split[i])
					((InterpSquare<V>*)(this->subTrees[i]))->updateCoeffs(d-1);
		}
	}
	
	/* virtual void merge_in(InterpSquare<V>* I)
	{
		
		for(int i=0; i<4; i++)
		{
			if(I->split[i])
			{
				if(split[i])
					(InterpSquare<V>*)(this->subTrees[i])->mergeIn((InterpSquare<V>*)I->subTrees[i]);
				else
				{
					delete(subTrees[i]);
					subTrees[i] = (InterpSquare<V>*)I->subTrees[i]
				}
			}
		}
		
	} */
};





template<class V>
class InterpGuardSquare: public InterpSquare<V>
{
	
public:
	/// Constructor
	InterpGuardSquare(TopMesh<V>* TM, const int d, const int ax, const int ay, V* v):
	InterpSquare<V>(TM, d, ax, ay, v)
	{ }
	
	/// Destructor
	virtual ~InterpGuardSquare() { }
	
	virtual void checkSubdivide(int d)
	{
		if(d)
		{
			for(int i=0; i<4; i++)
				if(this->split[i])
					this->subTrees[i]->checkSubdivide(d-1);
			return;
		}
		
		for(int dx=0; dx<2; dx++)
		{
			for(int dy=0; dy<2; dy++)
			{
				if(!dx && !dy) continue;
				V actualval = this->tm->evalAt(this->mydepth+1,2*this->absx+dx,2*this->absy+dy);
				setVal(1,dx,dy,actualval);
				V interpval = this->shallowInterpolate(1,dx,dy);
			}
		}		
	}
	
	virtual int doSplit(const int d)
	{
		int nsplit=0;
		
		if(!d)
		{
			for(int i = 0; i<4; i++)
			{
				if(this->split[i])
				{
					InterpGuardSquare* s = new InterpGuardSquare<V>(this->tm, (this->mydepth)+1, 2*(this->absx)+i%2, 2*(this->absy)+i/2, this->dtGetVal(1,i%2,i/2));
					this->subTrees[i]->forgetVal();
					delete(this->subTrees[i]);
					this->subTrees[i] = s;
					nsplit++;
				}
			}
			return nsplit;
		}
		
		for(int i = 0; i<4; i++)
			nsplit += this->subTrees[i]->doSplit(d-1);
		return nsplit;
	}
	
};


//TopMesh geometry is fundamentally toroidal; add non-splitting guard squares at margins to get a tube
template<class V>
class TopMesh
{
public:
	/// Constructor
	TopMesh(int NX, int NY, double X0, double Y0, double X1, double Y1, V (*ff)(double, double), InterpoBasis<V>* B, double relerr, double abserr):
	nx(NX), ny(NY), f(ff), b(B), xdivs((double*)malloc((NX+1)*sizeof(double))), ydivs((double*)malloc((NY+1)*sizeof(double))),
	squares((InterpSquare<V>**)malloc(NX*NY*sizeof(InterpSquare<V>*))), splitdepth(0), relerr2(relerr*relerr), abserr2(abserr*abserr)
	{
		for(int i=0; i<=nx; i++) xdivs[i] = X0 + (X1-X0)*double(i)/double(nx);
		for(int i=0; i<=ny; i++) ydivs[i] = Y0 + (Y1-Y0)*double(i)/double(ny);
		
		for(int i=0; i<nx; i++)
		{
			for(int j=0; j<ny; j++)
			{
				V* a = new V;
				*a = evalAt(0,i,j);
				squares[i + nx*j] = new InterpSquare<V>(this, 0, i, j, a);
			}
		}
		
		updateCoeffs(0);
		while(refine()) {}
	}
	
	TopMesh(int NX, int NY, double X0, double Y0, double X1, double Y1, V (*ff)(double, double), InterpoBasis<V>* B, double relerr, double abserr, bool NOSQUARES):
	nx(NX), ny(NY), f(ff), b(B), xdivs((double*)malloc((NX+1)*sizeof(double))), ydivs((double*)malloc((NY+1)*sizeof(double))),
	squares((InterpSquare<V>**)malloc(NX*NY*sizeof(InterpSquare<V>*))), splitdepth(0), relerr2(relerr*relerr), abserr2(abserr*abserr)
	{
		for(int i=0; i<=nx; i++) xdivs[i] = X0 + (X1-X0)*double(i)/double(nx);
		for(int i=0; i<=ny; i++) ydivs[i] = Y0 + (Y1-Y0)*double(i)/double(ny);
	}
	
	/// Destructor
	virtual ~TopMesh() { }
		
	virtual int refine()
	{
		checkSubdivide(splitdepth);
		int nsplit = doSplit(splitdepth);
		updateCoeffs(++splitdepth);
		printf("Mesh refined (%i new units)!\n",nsplit);
		return nsplit;
	}

	double absLocation(const int d, const int n, const bool xaxis) const
	{
		int c = n/(1<<d);
		int r = n%(1<<d);
		double* divs = xdivs;
		if(!xaxis) divs = ydivs;
		return divs[c]+(divs[c+1]-divs[c])*double(r)/double(1<<d);
	}
					  
	virtual void checkSubdivide(int d)
	{
		for(int i=0; i<nx*ny; i++)
			squares[i]->checkSubdivide(d);
	}
	
	virtual int doSplit(int d)
	{
		int nsplit = 0;
		for(int i=0; i<nx*ny; i++)
			nsplit += squares[i]->doSplit(d);
		return nsplit;
	}
	
	virtual double interpolate(double lx, double ly) const
	{
		return 0;
	}
	
	virtual void visualize(bool istop = true) const
	{
		if(istop) { vsr::startRecording(true); vsr::clearWindow(); }
		for(int i=0; i<nx*ny; i++)
			squares[i]->visualize(false);
		if(istop) vsr::stopRecording();
	}
	
	virtual void display(int v=0)
	{
		printf("TopMesh %ix%i\n",nx,ny);
		for(int i=0; i<nx*ny; i++) squares[i]->display(v);
	}
	
	virtual void updateCoeffs(int d)
	{
		for(int i=0; i<nx*ny; i++)
			squares[i]->updateCoeffs(d);
	}
	
	bool closeEnough(V actualval, V interpval) const
	{
		double m = (actualval-interpval).mag2();
		if(m < abserr2) return true;
		if(m/actualval.mag2() < relerr2) return true;
		return false;
	}
	
	V* getVal(const int d, const int lx, const int ly)
	{
		return squares[mainSquare(d,lx,ly)]->getVal(d,(lx+(1<<d+5))%(1<<d),(ly+(1<<d+5))%(1<<d));
	}
	
	V evalAt(const int d, const int lx, const int ly) const
	{
		return f(absLocation(d,lx,true),absLocation(d,ly,false));
	}
	
	void markSplitAffected(const int d, const int lx, const int ly)
	{
		for(int i = 0; i < b->nsupport; i++)
			squares[mainSquare(d,lx-b->supportx[i],ly-b->supporty[i])]->markForSplit(d,(lx-b->supportx[i]+(1<<d+5))%(1<<d),(ly-b->supporty[i]+(1<<d+5))%(1<<d));
	}
	
	int mainSquare(const int d, const int lx, const int ly) const
	{
		int s = 1<<d;
		int rx = ((lx+9*nx*s)%(nx*s))/s;
		int ry = ((ly+9*ny*s)%(ny*s))/s;
		return rx + nx*ry;
	}
	
	InterpoBasis<V>* b;
	
protected:
	int nx,ny;
	int splitdepth;
	double* xdivs;
	double* ydivs;
	double relerr2;
	double abserr2;
	InterpSquare<V>** squares;
	V (*f)(double, double);
};


template<class V>
class TubeMesh: public TopMesh<V>
{
public:
	/// Constructor
	TubeMesh(int NX, int NY, double Y0, double Y1, V (*ff)(double, double), InterpoBasis<V>* B, double relerr, double abserr):
	TopMesh<V>(NX, NY+2, -M_PI, Y0-(Y1-Y0)/double(NY), M_PI, Y1+(Y1-Y0)/double(NY), ff, B, relerr, abserr, true)
	{
		for(int i=0; i<this->nx; i++)
		{
			for(int j=0; j<this->ny; j++)
			{
				V* a = new V;
				*a = this->evalAt(0,i,j);
				if(j == 0 || j == this->ny-1)
					this->squares[i + this->nx*j] = new InterpGuardSquare<V>(this, 0, i, j, a);
				else
					this->squares[i + this->nx*j] = new InterpSquare<V>(this, 0, i, j, a);
			}
		}	
		
		this->updateCoeffs(0);
		while(this->refine()) { }
	}
	
	/// Destructor
	virtual ~TubeMesh() { }
	
	
	virtual void visualize(bool istop = true) const
	{
		if(istop) { vsr::startRecording(true); vsr::clearWindow(); }
		for(int i=0; i<this->nx*this->ny; i++)
			this->squares[i]->visualize_cyl();
		if(istop) vsr::stopRecording();
	}
};


#endif
