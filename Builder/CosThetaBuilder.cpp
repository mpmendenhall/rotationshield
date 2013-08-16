#include "CosThetaBuilder.hh"
#include "strutils.hh"

mdouble ShiftPositioner::angle(unsigned int i, unsigned int ncoils) const {
	i = i%ncoils;
	mdouble x0 = (0.5+i)/(mdouble)ncoils;
	mdouble x = x0;
	for(unsigned int k=0; k<shift.size(); k++) x += shift[k]*sin(PI*(mdouble)(k+1)*x0);
	mdouble y = sqrt(1-x*x);
	return atan2(y,x);
}

Stringmap ShiftPositioner::getInfo() const {
	Stringmap m;
	m.insert("class","ShiftPositioner");
	return m;
}

mdouble shiftPositioner(unsigned int i, unsigned int ncoils, void* params) {
	i = i%ncoils;
	mdouble x0 = (0.5+i)/(mdouble)ncoils;
	mdouble x = x0;
	VarVec<mdouble>* shift = (VarVec<mdouble>*)params;
	if(shift)
		for(unsigned int k=0; k<shift->size(); k++) x += (*shift)[k]*sin(PI*(mdouble)(k+1)*x0);
	mdouble y = sqrt(1-x*x);
	return atan2(y,x);
}

mdouble alarconKPositioner(unsigned int i, unsigned int ncoils, void* params) {	
	i = i%ncoils;
	mdouble k = *(mdouble*)params;
	mdouble x0 = (0.5+i)/(mdouble)ncoils;
	return atan2(sqrt(1/(x0*x0)-1),(1-k)*(1-k));
}


//--------------------------------------------------------------------

Stringmap VecTrans::getInfo() const {
	Stringmap m;
	m.insert("class","VecTrans");
	for(unsigned int i=0; i<3; i++) {
		m.insert(std::string("v1_")+itos(i),trans1[i]);
		m.insert(std::string("v2_")+itos(i),trans2[i]);
	}
	return m;
}

//--------------------------------------------------------------------

void CosThetaBuilder::writeInfo(QFile& qOut) const {
	if(AP) qOut.insert("coilPositioner",AP->getInfo());
	if(ET) qOut.insert("coilTrans",ET->getInfo());
	Stringmap m;
	m.insert("length",length);
	m.insert("radius",radius);
	m.insert("ncoils",itos(ncoils));
	qOut.insert("cosThetaCoil",m);
}

void CosThetaBuilder::buildCoil(MixedSource& M) {
	buildEndpoints();
	buildSides(M);
	buildLineCaps(M);
}

vec3 CosThetaBuilder::getEndp(unsigned int n, bool xside, bool yside, bool zside) {
	assert(n<ncoils);
	return endp[8*n + zside + 2*yside + 4*xside];
}

void CosThetaBuilder::buildEndpoints() {
	assert(AP);
	mdouble x,y,phi;
	mdouble xmul,ymul,zmul;
	endp = std::vector<vec3>(8*ncoils);
	for(unsigned int i = 0; i<ncoils; i++) {
		for(int a1 = 0; a1<2; a1++) { // z
			for(int a2=0; a2<2; a2++) { // y
				for(int a3=0; a3<2; a3++) { // x
					xmul = (a3)?-1:1;
					ymul = (a2)?-1:1;
					zmul = (a1)?-1:1;
					phi = AP->angle(i+a3*ncoils+a2*2*ncoils+a1*4*ncoils,ncoils);
					x = radius*cos(phi);
					y = radius*sin(phi);
					endp[8*i + a1 + 2*a2 + 4*a3] = vec3(x*xmul,y*ymul,zmul*length*0.5) + (ET?ET->trans(ncoils,i,a1,a2,a3):vec3());
				}
			}
		}
	}
}

void CosThetaBuilder::buildSides(MixedSource& M) {
	for(unsigned int i=0; i<ncoils; i++)
		for(unsigned int xside = 0; xside < 2; xside++)
			for(unsigned int yside = 0; yside < 2; yside++)
				M.addsource(new LineSource(getEndp(i,xside,yside,!yside),getEndp(i,xside,yside,yside),1.0/mdouble(ncoils)));
}

void CosThetaBuilder::buildArcCaps(MixedSource& M,unsigned int nseg) {
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int zside = 0; zside < 2; zside++) {
			
			for(unsigned int yside = 0; yside < 2; yside++) {
				for(unsigned int i=0; i<ncoils-1; i++) {
					if(yside == zside)
						M.arc(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)/mdouble(ncoils),nseg);
					else
						M.arc(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)/mdouble(ncoils),nseg);
				}
			}
			
			M.arc(getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),1.0,nseg);
		}
	}
	
}

void CosThetaBuilder::buildLineCaps(MixedSource& M) {
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int zside = 0; zside < 2; zside++) {
			
			for(unsigned int yside = 0; yside < 2; yside++) {
				for(unsigned int i=0; i<ncoils-1; i++) {
					if(yside == zside)
						M.addsource( new LineSource(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)/mdouble(ncoils)) );
					else
						M.addsource( new LineSource(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)/mdouble(ncoils)) );
				}
			}
			
			M.addsource( new LineSource( getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),1.0) );
		}
	}
	
}

void CosThetaBuilder::buildStraightCaps(MixedSource& M) {
	for(unsigned int xside = 0; xside < 2; xside++)
		for(unsigned int zside = 0; zside < 2; zside++)
			for(unsigned int i=0; i<ncoils; i++)
				M.addsource(new LineSource(getEndp(i,xside,zside,zside),getEndp(i,xside,!zside,zside),1.0/mdouble(ncoils)));
}

void CosThetaBuilder::buildMixedCaps(MixedSource& M, float rinner) {
	
	std::vector<vec3> innercircle[2][2][2];
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int zside = 0; zside < 2; zside++) {
			for(unsigned int i=0; i<ncoils; i++) {
				vec3 v1 = getEndp(i,xside,zside,zside);
				vec3 v2 = getEndp(i,xside,!zside,zside);
				if(fabs(v1[0]) >= rinner*radius) {
					M.addsource(new LineSource(v1,v2,1.0/mdouble(ncoils)));
				} else {
					float l = (1-sqrt(rinner*radius*rinner*radius-v1[0]*v1[0])/radius)*0.5;
					M.addsource(new LineSource(v1,v1*(1-l) + v2*l,1.0/mdouble(ncoils)));
					M.addsource(new LineSource(v1*l + v2*(1-l),v2,1.0/mdouble(ncoils)));
					innercircle[xside][0][zside].push_back(v1*(1-l) + v2*l);
					innercircle[xside][1][zside].push_back(v1*l + v2*(1-l));
				}
			}
		}
	}
	
	unsigned int imax = innercircle[0][0][0].size();
	vec3 dz[2];
	dz[0] = vec3(0,0,0.20*length);
	dz[1] = -dz[0];
	
	if(imax) {
		for(unsigned int xside = 0; xside < 2; xside++) {
			for(unsigned int zside = 0; zside < 2; zside++) {
				for(unsigned int yside = 0; yside < 2; yside++) {
					
					for(unsigned int i=0; i<imax-1; i++) {
						if(yside)
							M.arc(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i+1]+dz[zside],(i+1)/mdouble(ncoils),32);
						else
							M.arc(innercircle[xside][yside][zside][i+1]+dz[zside],innercircle[xside][yside][zside][i]+dz[zside],(i+1)/mdouble(ncoils),32);
					}
					
					for(unsigned int i=0; i<imax; i++) {
						if(!yside)
							M.addsource(new LineSource(innercircle[xside][yside][zside][i],innercircle[xside][yside][zside][i]+dz[zside],1/mdouble(ncoils)));
						else
							M.addsource(new LineSource(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i],1/mdouble(ncoils)));
					}
				}
				M.arc(innercircle[xside][true][zside][imax-1]+dz[zside],innercircle[xside][false][zside][imax-1]+dz[zside],mdouble(imax)/mdouble(ncoils),32);
			}
		}
	}
	
}


