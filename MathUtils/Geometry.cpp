#include "Geometry.hh"

void Line::visualize(bool top) const { 
	if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); } 
	vsr::Visr::W->line(s,e); 
	if(top) vsr::Visr::W->stopRecording();
}

void Line::visualizeDirected(float j) const { 
	vsr::Visr::W->startLines();
	for(float i=0; i<8; i++) {
		float p = i/7.0;
		float q = 1.0-p;
		if(p<0.5)
			vsr::Visr::W->setColor(2*p,0,1.0,1.0);
		else
			vsr::Visr::W->setColor(1.0,0,2*q,1.0);
		if(j>0)
			vsr::Visr::W->vertex(s*q+e*p);
		else
			vsr::Visr::W->vertex(s*p+e*q);
	}
	vsr::Visr::W->endLines();
}

Plane::Plane(annulusSpec a) {
	vec2 c = (a.start+a.end)*0.5;
	vec2 d = a.end-a.start;
	o = zrotate(vec3(0,c[1],c[0]),a.theta0);
	dx = zrotate(vec3(2*c[1]*tan(a.dTheta*0.5),0,0),a.theta0);
	dz = zrotate(vec3(0,d[1],d[0]),a.theta0);
	wx = dx.mag();
	wz = dz.mag();
	area = wx*wz;
	sn = cross(dx,dz)/area;
}

const mat3 Plane::projectionMatrix() const {
	mat3 M;
	for(unsigned int i=0; i<3; i++) {
		M(0,i) = dx[i]/wx;
		M(1,i) = dz[i]/wz;
		M(2,i) = sn[i];
	}
	return M;
}


void Plane::visualize(bool top) const { 
	if(top) { vsr::Visr::W->startRecording(); vsr::Visr::W->clearWindow(); } 
	vsr::Visr::W->plane(o,dx,dz);
	if(top) vsr::Visr::W->stopRecording();
}

void Plane::visualizeCoordinates(float scale) const {
	float l = min(wx,wz)*0.5*scale;
	vsr::Visr::W->line(o, o+dx*l/wx);
	vsr::Visr::W->line(o, o+dz*l/wz);
	vsr::Visr::W->line(o, o+sn*l);
}

void Plane::visualizeVector(vec3 v) const {
	v = v.normalized()*(0.5*min(wx,wz));
	vsr::Visr::W->line(o,o+dx*v[0]/wx+dz*v[1]/wz+sn*v[2]);
}

Plane Plane::readFromFile(std::istream& s) { 
	vec3 origin = vec3::readBinary(s);
	vec3 dxx = vec3::readBinary(s);
	vec3 dzz = vec3::readBinary(s);
	return Plane(origin,dxx,dzz);
}


vec3 Plane::zrotate(const vec3& v, mdouble th) {
	mdouble cs = cos(th);
	mdouble sn = sin(th);
	return vec3( cs*v[0]-sn*v[1], sn*v[0]+cs*v[1], v[2] );
}

std::vector<Plane> Plane::subdivide(unsigned int nX, unsigned int nZ) const {
	assert(nX>0 && nZ>0);
	std::vector<Plane> v = std::vector<Plane>();
	for(unsigned int x=0; x<nX; x++)
		for(unsigned int z=0; z<nZ; z++)
			v.push_back( Plane(o - dx*(0.5 - (mdouble(x)+0.5)/nX) - dz*(0.5 - (mdouble(z)+0.5)/nZ), dx/nX, dz/nZ) );
	return v;
}