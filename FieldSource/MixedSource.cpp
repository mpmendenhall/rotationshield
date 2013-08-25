#include "MixedSource.hh"
#include "ScanRange.hh"

void MixedSource::loadSourcesFile(FILE* f, mdouble scale = 1.0) {
	float x1,y1,z1,x2,y2,z2,j;
	while(fscanf(f,"LINE %f %f %f %f %f %f %f ",&x1,&y1,&z1,&x2,&y2,&z2,&j) == 7)
		addsource(new LineSource(vec3(x1,y1,z1),vec3(x2,y2,z2),j*scale));
}

vec3 MixedSource::fieldAt(const vec3& v) const {
	vec3 b(0,0,0);
	for(unsigned int i=0; i<sources.size(); i++) b += sources[i]->fieldAt(v);
	return b;	
}

vec3 MixedSource::fieldOverLine(Line l) const {
	vec3 b = vec3();
	for(unsigned int i=0; i<sources.size(); i++) b += sources[i]->fieldOverLine(l);
	return b;	
}

vec3 MixedSource::fieldOverPlane(Plane p) const {
	vec3 b = vec3();
	for(unsigned int i=0; i<sources.size(); i++) b += sources[i]->fieldOverPlane(p);
	return b;	
}

void MixedSource::arc(vec3 start, vec3 end, mdouble j, int nsegs) {
	mdouble r = sqrt(start[0]*start[0]+start[1]*start[1]);
	mdouble z = start[2];
	mdouble th0 = atan2(start[1],start[0]);
	mdouble th1 = atan2(end[1],end[0]);
	if(th1-th0>M_PI) th0 += 2*M_PI;
	if(th0-th1>M_PI) th1 += 2*M_PI;
	if(th0 > th1)
	{
		mdouble tmp = th1; th1 = th0; th0 = tmp;
		j = -j;
	}
	int nsteps = int((th1-th0)/(2.0*M_PI/mdouble(nsegs))+1);
	mdouble th = th0;
	vec3 v1,v2; v1[2] = v2[2] = z;
	v1[0] = r*cos(th); v1[1] = r*sin(th);
	for(int i=1; i<nsteps+1; i++)
	{
		th = th0+mdouble(i)/mdouble(nsteps)*(th1-th0);
		v2[0] = r*cos(th); v2[1] = r*sin(th);
		addsource(new LineSource(v1,v2,j));
		v1 = v2;
	}
}

void MixedSource::loop(mdouble z, mdouble r, mdouble j, int nsides) {
	vec3 v0,v1;
	v0[2] = v1[2] = z;
	v0[0] = r; v0[1] = 0;
	ScanRange sr(0,2*M_PI,nsides);
	mdouble th = sr.next();
	v0[0] = r*cos(th); v0[1] = r*sin(th);
	for(th = sr.next(); sr.goOn(); th=sr.next())
	{
		v1[0] = r*cos(th); v1[1] = r*sin(th);
		addsource(new LineSource(v0,v1,j));
		v0 = v1;
	}
}
