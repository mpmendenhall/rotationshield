#include "PlanarElement.hh"
#include <algorithm>

void PlanarElement::_visualize() const {

	mdouble logmax = 3.0;
	float smag = state.mag2();
	if(smag)
		smag = std::max(0.,std::min(1.0,0.1*(log(smag)+10-logmax)));
	
	vec3 svec = vec3();
	for(unsigned int i=0; i<std::min((unsigned int)3, state.size()); i++)
		svec[i] = state[i];
	
	vec3 hsv = vec3( atan2(svec[0],svec[1]), 1.0, 1.0 );
	vec3 rgb = hsv2rgb(hsv);
	vsr::setColor(rgb[0], rgb[1] , rgb[2], 1.0*smag);
	p._visualize();
	vsr::setColor(1, .7, .7, 1.0);
	p.visualizeCoordinates(0.2);
	vsr::setColor(.7, .7, 1, 1.0);
	p.visualizeCoordinates(-0.2);
	vsr::setColor(1,1,1,1);
	
	p.visualizeVector(svec);
}
