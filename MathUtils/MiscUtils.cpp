#include "MiscUtils.hh"
#include <cstdio>

mdouble randunif(mdouble a, mdouble b) { 
	return a + (b-a)*mdouble(rand())/mdouble(RAND_MAX);
}

float normalizeAngle(float a, float theta0) {
	if(a<theta0) a += 2*PI*int(1+(theta0-a)/(2*PI));
	if(a>=theta0+2*PI) a -= 2*PI*int((a-theta0)/(2*PI));
	return a;
}
