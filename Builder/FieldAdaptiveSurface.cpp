#include "FieldAdaptiveSurface.hh"

FieldAdaptiveSurface::FieldAdaptiveSurface(const DVFunc1<2,mdouble>& f): F(f) {
	unsigned int npts = 100;
	Ilz.setupDataGrid(&npts,&npts+1);
	Ilz.getInterpolator()->scale = float(npts+1)/float(npts);
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
	const unsigned int ngridpts = Ilz.n_pts();
	double a[ngridpts];
	for(unsigned int i=0; i<ngridpts; i++)
		a[i] = float(i)/float(ngridpts-1);
	Ilz.setData(a);
}

void FieldAdaptiveSurface::optimizeSpacing(const FieldEstimator2D& fes, double pfixed) {
		
	const unsigned int ngridpts = Ilz.n_pts();
	
	// cumulative parameter across length
	double fstr[ngridpts];
	fstr[0]=0;
	for(unsigned int i=1; i<ngridpts; i++) {
		mdouble l = float(i)/(ngridpts-1.);
		vec2 zr = F(l);
		vec2 dv = F.deriv(l);
		// bunch up at high field derivatives along profile curve
		fstr[i] = fstr[i-1] + fabs(fes.derivAt(zr,dv).dot(dv));
	}
	
	// add "slope" to cumulative curve for fixed partitions
	pfixed = pfixed<=1e-3 ? 1e-3 : (pfixed>=0.999 ? 0.999:pfixed);
	mdouble slope = (fstr[ngridpts-1]+1e-8) * pfixed/(1.-pfixed) / float(ngridpts);
	printf("Optimized surface z spacing: to accumulated x=%g, adding slope %g\n",fstr[ngridpts-1],slope*ngridpts);
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] += i*slope;
	
	// convert to cumulative curve; set interpolator
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] /= fstr[ngridpts-1];
	Ilz.setData(fstr);
}


