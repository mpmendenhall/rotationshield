#include "ShieldBuilder.hh"
#include <math.h>

vec2 FieldEstimator2D::estimateAt(const vec2& v) const {
	vec2 est = vec2();
	vec2 r;
	mdouble r0;
	for(unsigned int i=0; i<sources.size(); i++) {
		r = v-sources[i];
		r0 = r.mag2();
		est[0] += r[1]*(-currents[i]/r0);
		est[1] += r[0]*(currents[i]/r0);
	}
	return est;
}

vec2 FieldEstimator2D::derivAt(const vec2& v, const vec2& dx) const {
	return (estimateAt(v+dx*0.5)-estimateAt(v-dx*0.5))/dx.mag();
}

vec2 FieldEstimator2Dfrom3D::estimateAt(const vec2& v) const {
	vec3 F = myFS->fieldAt(vec3(v[1],0,v[0]));
	return vec2(F[2],F[0]);
}

//------------------------------------------------------

void ShieldBuilder::OptCone(unsigned int nZ0, unsigned int nZ1, vec2 s, vec2 e, PlanarElement* base, FieldEstimator2D* fes ) {
	
	printf("Optimizing shield grid with %i fixed segments, %i varying segments\n",nZ0,nZ1);
	if(fes) printf("Field estimator provided.\n");
	
	const unsigned int ngridpts = 2001;
	const unsigned int nZ = nZ0+nZ1;
	
	//cumulative field strength across length
	mdouble fstr[ngridpts];
	fstr[0]=0;
	vec2 dv = (e-s)*1./ngridpts;
	if(fes && nZ1) {
		for(unsigned int i=1; i<ngridpts; i++) {
			float l = (float(i)-0.5)/(ngridpts-1.0);
			//fstr[i] = fstr[i-1] + fes->estimateAt(s*(1-l)+e*l).mag(); 		// bunch up at high fields
			fstr[i] = fstr[i-1] + fabs(fes->derivAt(s*(1-l)+e*l,dv).dot(dv));	// bunch up at high field derivatives along sampling
		}
	} else {
		for(unsigned int i=0; i<ngridpts; i++)
			fstr[i] = i;
	}
	
	printf("Average field %g\n",fstr[ngridpts-1]/(ngridpts-1));
	
	mdouble slope = 0;
	if(nZ1)
		slope = fstr[ngridpts-1]*nZ0/((ngridpts-1.0)*nZ1); // add a constant slope for fixed partitions
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] += i*slope;
	
	for(unsigned int i=0; i<ngridpts; i++)
		fstr[i] *= nZ/fstr[ngridpts-1];
	
	//interpolation to determine dividing lines, relative coordinates 0 to 1
	mdouble* ls = new mdouble[nZ+1];
	unsigned int n=1;
	int i=1;
	ls[0]=0;
	while(n<nZ) {
		if(fstr[i] >= n) {
			ls[n] = (mdouble(i)-(fstr[i]-n)/(fstr[i]-fstr[i-1]))/(ngridpts-1.0); // interpolated partition point
			n++;
			continue;
		}
		i++;	
	}
	ls[nZ] = 1.0;
	
	annulusSpec a;
	a.theta0 = 0;
	a.dTheta = 2*M_PI/mdouble(nTheta);
	
	base->retain();
	for(unsigned int z=0; z<nZ; z++) {
		a.start = s*(1-ls[z])+e*ls[z];
		a.end = s*(1-ls[z+1])+e*ls[z+1];
		append(new ShieldSegment(nTheta, base->reference(a)));
	}
	base->release();
	
	delete(ls);
}
