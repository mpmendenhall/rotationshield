#include "FieldAdaptiveSurface.hh"

FieldAdaptiveSurface::FieldAdaptiveSurface(const DVFunc1<2,mdouble>& f): F(f) {
	unsigned int npts = 100;
	Ilz.setupDataGrid(&npts,&npts+1);
	Ilz.getInterpolator()->scale = float(npts)/float(npts-3);
	Ilz.getInterpolator()->offset = -1./float(npts-3);
	setConstantSpacing();
}

vec2 FieldAdaptiveSurface::operator()(mdouble x) const {
	double l = x;
	return F(Ilz.eval(&l));
}

mdouble FieldAdaptiveSurface::l_dist_deriv(mdouble l) const {
	const double h = 1e-3;
	double l0 = l-h;
	double l1 = l+h;
	return (Ilz.eval(&l1)-Ilz.eval(&l0))/(2*h);
}
	
vec2 FieldAdaptiveSurface::deriv(mdouble x) const {
	double l = x;
	return F.deriv(Ilz.eval(&l)) * l_dist_deriv(x);
}

void FieldAdaptiveSurface::setConstantSpacing() {
	const unsigned int npts = Ilz.n_pts();
	double a[npts];
	for(unsigned int i=0; i<npts; i++)
		a[i] = (double(i)-1.)/double(npts-3);
	Ilz.setData(a);
}

void FieldAdaptiveSurface::symmetry_test() const {
	for(int i=0; i<21; i++) {
		double l = float(i)/20.;
		double v =Ilz.eval(&l);
		printf("l=%g:\tv=%g\t%g\n",l,v,1-v);
	}
}

void FieldAdaptiveSurface::optimizeSpacing(const FieldEstimator2D& fes, double pfixed) {
	
	const unsigned int npts = Ilz.n_pts();
	
	// cumulative parameter across length
	double fstr[npts];
	fstr[1]=0;
	for(unsigned int i=2; i <= npts-2; i++) {
		mdouble l = float(i-1.5)/(npts-3.);
		vec2 zr = F(l);
		vec2 dv = F.deriv(l);
		// bunch up at high field derivatives along profile curve
		fstr[i] = fstr[i-1] + fabs(fes.derivAt(zr,dv).dot(dv));
	}
	
	// add "slope" to cumulative curve for fixed partitions
	pfixed = pfixed<=1e-3 ? 1e-3 : (pfixed>=0.999 ? 0.999:pfixed);
	mdouble slope = (fstr[npts-2]+1e-8) * pfixed/(1.-pfixed) / float(npts);
	printf("Optimized surface z spacing: to accumulated x=%g, adding slope %g\n",fstr[npts-2],slope*npts);
	for(unsigned int i=2; i <= npts-2; i++)
		fstr[i] += (i-1)*slope;
	
	// convert to cumulative curve;
	for(unsigned int i=2; i<=npts-2; i++)
		fstr[i] /= fstr[npts-2];
	// special endpoints for linear approach
	fstr[0] = -fstr[2];
	fstr[npts-1] = 2-fstr[npts-3];
	// set interpolator
	Ilz.setData(fstr);
}


