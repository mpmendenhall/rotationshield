#include "FieldEstimator2D.hh"

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
