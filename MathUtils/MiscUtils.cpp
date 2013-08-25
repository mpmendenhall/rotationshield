#include "MiscUtils.hh"
#include <cstdio>
#include <math.h>

mdouble randunif(mdouble a, mdouble b) { 
	return a + (b-a)*mdouble(rand())/mdouble(RAND_MAX);
}

float normalizeAngle(float a, float theta0) {
	if(a<theta0) a += 2*M_PI*int(1+(theta0-a)/(2*M_PI));
	if(a>=theta0+2*M_PI) a -= 2*M_PI*int((a-theta0)/(2*M_PI));
	return a;
}
