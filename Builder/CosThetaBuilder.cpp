/* 
 * CosThetaBuilder.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "CosThetaBuilder.hh"
#include "strutils.hh"
#include <math.h>

double ShiftPositioner::angle(unsigned int i, unsigned int ncoils) const {
	i = i%ncoils;
	double x0 = (0.5+i)/(double)ncoils;
	double x = x0;
	for(unsigned int k=0; k<shift.size(); k++) x += shift[k]*sin(M_PI*(double)(k+1)*x0);
	double y = sqrt(1-x*x);
	return atan2(y,x);
}

Stringmap ShiftPositioner::getInfo() const {
	Stringmap m;
	m.insert("class","ShiftPositioner");
	if(shift.size())
		m.insert("shift",vtos(varvec2doublevec<double>(shift)));
	return m;
}

//--------------------------------------------------------------------

double AlarconKPositioner::angle(unsigned int i, unsigned int ncoils) const {
	i = i%ncoils;
	double x0 = (0.5+i)/(double)ncoils;
	return atan2(sqrt(1/(x0*x0)-1),(1-k)*(1-k));
}

Stringmap AlarconKPositioner::getInfo() const {
	Stringmap m;
	m.insert("class","AlarconKPositioner");
	m.insert("k",k);
	return m;
}

//--------------------------------------------------------------------

Stringmap VecTrans::getInfo() const {
	Stringmap m;
	m.insert("class","VecTrans");
	m.insert(std::string("t1"),vtos(vec2doublevec<3,double>(trans1)));
	m.insert(std::string("t2"),vtos(vec2doublevec<3,double>(trans2)));
	return m;
}

//--------------------------------------------------------------------

void mi_setGeometry(StreamInteractor* S) {
	float j = S->popFloat();
	float l = S->popFloat();
	float r = S->popFloat();
	float n = S->popInt();
	
	CosThetaBuilder* CT = dynamic_cast<CosThetaBuilder*>(S);
	CT->ncoils = n;
	CT->radius = r;
	CT->length = l;
	CT->j_total = j;
}

void mi_setDistortion(StreamInteractor* S) {
	float a = S->popFloat();
	int n = S->popInt();
	
	CosThetaBuilder* CT = dynamic_cast<CosThetaBuilder*>(S);
	assert(CT->AP);
	ShiftPositioner* SP = dynamic_cast<ShiftPositioner*>(CT->AP);
	if(!SP) {
		printf("Changing to ShiftPositioner type distortion\n");
		delete CT->AP;
		SP = new ShiftPositioner();
		CT->AP = SP;
	}
	if(n<=0) SP->shift = VarVec<double>(0);
	else {
		while(SP->shift.size()<(unsigned int)n) SP->shift.push_back(0);
		SP->shift[n-1] = a;
	}
}

void mi_setDistortionK(StreamInteractor* S) {
	float k = S->popFloat();
	
	CosThetaBuilder* CT = dynamic_cast<CosThetaBuilder*>(S);
	assert(CT->AP);
	delete CT->AP;
	CT->AP = new AlarconKPositioner(k);
}


void mi_setEndcaps(StreamInteractor* S) {
	std::string ec[2];
	ec[1] = S->popString();
	ec[0] = S->popString();
	CosThetaBuilder* CT = dynamic_cast<CosThetaBuilder*>(S);
	for(unsigned int i=0; i<2; i++) {
		if(ec[i]=="arc") CT->myCap[i] = CosThetaBuilder::CAP_ARC;
		if(ec[i]=="line") CT->myCap[i] = CosThetaBuilder::CAP_LINE;
		if(ec[i]=="none") CT->myCap[i] = CosThetaBuilder::CAP_NONE;
	}
}

void mi_setTranslation(StreamInteractor* S) {
	float z = S->popFloat();
	float y = S->popFloat();
	float x = S->popFloat();
	CosThetaBuilder* CT = dynamic_cast<CosThetaBuilder*>(S);
	CT->ET = new VecTrans(vec3(x,y,z));
}

//--------------------------------------------------------------------

CosThetaBuilder::CosThetaBuilder(unsigned int n, double r, double l, double j, AnglePositioner* ap, EndTranslator* et):
ncoils(n), radius(r), length(l), j_total(j), AP(ap), ET(et), nArc(100),
setGeometry("Geometry", &mi_setGeometry, this),
setDistortion("Distortion", &mi_setDistortion, this),
setDistortionK("Distortion K", &mi_setDistortionK, this),
selectEndcap("end wire shape"),
setEndcaps("End wires", &mi_setEndcaps, this),
setTranslation("Position offset", &mi_setTranslation, this),
OMcoil("Cos Theta Coil") {

	myCap[0] = myCap[1] = CAP_ARC;

	setGeometry.addArg("half-N loops","15");
	setGeometry.addArg("radius","0.5");
	setGeometry.addArg("length","1");
	setGeometry.addArg("current","1");
	//
	setDistortion.addArg("n","1");
	setDistortion.addArg("a_n","-0.001");
	//
	setDistortionK.addArg("k","0.005");
	//
	selectEndcap.addChoice("smooth arc","arc");
	selectEndcap.addChoice("straight line","line");
	selectEndcap.addChoice("no end wires","none");
	selectEndcap.setDefault("arc");
	//
	setEndcaps.addArg(&selectEndcap,"-z");
	setEndcaps.addArg(&selectEndcap,"+z");
	//
	setTranslation.addArg("x","0");
	setTranslation.addArg("y","0");
	setTranslation.addArg("z","0");
	//
	OMcoil.addChoice(&setGeometry,"geom");
	OMcoil.addChoice(&setDistortion,"dist");
	OMcoil.addChoice(&setDistortionK,"kdist");
	OMcoil.addChoice(&setEndcaps,"ends");
	OMcoil.addChoice(&setTranslation,"off");
	
}


void CosThetaBuilder::writeInfo(QFile& qOut) const {
	if(AP) qOut.insert("coilPositioner",AP->getInfo());
	if(ET) qOut.insert("coilTrans",ET->getInfo());
	Stringmap m;
	m.insert("length",length);
	m.insert("radius",radius);
	m.insert("ncoils",itos(ncoils));
	m.insert("current",j_total);
	qOut.insert("cosThetaCoil",m);
}

void CosThetaBuilder::buildCoil(MixedSource& M) {
	buildEndpoints();
	buildSides(M);
	for(unsigned int zside=0; zside<2; zside++) {
		if(myCap[zside]==CAP_LINE)
			buildLineCap(M,zside);
		else if(myCap[zside]==CAP_STRAIGHT)
			buildStraightCap(M,zside);
		else if(myCap[zside]==CAP_ARC)
			buildArcCap(M,zside,nArc);
		else if(myCap[zside]==CAP_NONE) {}
		else assert(false);
	}
}

vec3 CosThetaBuilder::getEndp(unsigned int n, bool xside, bool yside, bool zside) {
	assert(n<ncoils);
	return endp[8*n + zside + 2*yside + 4*xside];
}

void CosThetaBuilder::buildEndpoints() {
	assert(AP);
	double x,y,phi;
	double xmul,ymul,zmul;
	endp = std::vector<vec3>(8*ncoils);
	for(unsigned int i = 0; i<ncoils; i++) {
		for(int a1 = 0; a1<2; a1++) { // z
			for(int a2=0; a2<2; a2++) { // y
				for(int a3=0; a3<2; a3++) { // x
					xmul = (a3)?-1:1;
					ymul = (a2)?-1:1;
					zmul = (a1)?-1:1;
					phi = AP->angle(i + a3*ncoils + a2*2*ncoils + a1*4*ncoils, ncoils);
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
				M.addsource(new LineSource(getEndp(i,xside,yside,!yside),getEndp(i,xside,yside,yside), j_total/double(ncoils)));
}

void CosThetaBuilder::buildArcCap(MixedSource& M, unsigned int zside, unsigned int nseg) {
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int yside = 0; yside < 2; yside++) {
			for(unsigned int i=0; i<ncoils-1; i++) {
				if(yside == zside)
					M.arc(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)*j_total/double(ncoils),nseg);
				else
					M.arc(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)*j_total/double(ncoils),nseg);
			}
		}
		M.arc(getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),j_total,nseg);
	}
}

void CosThetaBuilder::buildLineCap(MixedSource& M, unsigned int zside) {
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int yside = 0; yside < 2; yside++) {
			for(unsigned int i=0; i<ncoils-1; i++) {
				if(yside == zside)
					M.addsource( new LineSource(getEndp(i,xside,yside,zside),getEndp(i+1,xside,yside,zside),(i+1)*j_total/double(ncoils)) );
				else
					M.addsource( new LineSource(getEndp(i+1,xside,yside,zside),getEndp(i,xside,yside,zside),(i+1)*j_total/double(ncoils)) );
			}
		}
		M.addsource( new LineSource( getEndp(ncoils-1,xside,zside,zside),getEndp(ncoils-1,xside,!zside,zside),j_total) );
	}
}

void CosThetaBuilder::buildStraightCap(MixedSource& M, unsigned int zside) {
	for(unsigned int xside = 0; xside < 2; xside++)
		for(unsigned int i=0; i<ncoils; i++)
			M.addsource(new LineSource(getEndp(i,xside,zside,zside),getEndp(i,xside,!zside,zside),j_total/double(ncoils)));
}

void CosThetaBuilder::buildMixedCaps(MixedSource& M, float rinner) {
	
	std::vector<vec3> innercircle[2][2][2];
	
	for(unsigned int xside = 0; xside < 2; xside++) {
		for(unsigned int zside = 0; zside < 2; zside++) {
			for(unsigned int i=0; i<ncoils; i++) {
				vec3 v1 = getEndp(i,xside,zside,zside);
				vec3 v2 = getEndp(i,xside,!zside,zside);
				if(fabs(v1[0]) >= rinner*radius) {
					M.addsource(new LineSource(v1,v2,j_total/double(ncoils)));
				} else {
					float l = (1-sqrt(rinner*radius*rinner*radius-v1[0]*v1[0])/radius)*0.5;
					M.addsource(new LineSource(v1,v1*(1-l) + v2*l,j_total/double(ncoils)));
					M.addsource(new LineSource(v1*l + v2*(1-l),v2,j_total/double(ncoils)));
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
							M.arc(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i+1]+dz[zside], (i+1)*j_total/double(ncoils),32);
						else
							M.arc(innercircle[xside][yside][zside][i+1]+dz[zside],innercircle[xside][yside][zside][i]+dz[zside], (i+1)*j_total/double(ncoils),32);
					}
					
					for(unsigned int i=0; i<imax; i++) {
						if(!yside)
							M.addsource(new LineSource(innercircle[xside][yside][zside][i],innercircle[xside][yside][zside][i]+dz[zside],j_total/double(ncoils)));
						else
							M.addsource(new LineSource(innercircle[xside][yside][zside][i]+dz[zside],innercircle[xside][yside][zside][i],j_total/double(ncoils)));
					}
				}
				M.arc(innercircle[xside][true][zside][imax-1]+dz[zside],innercircle[xside][false][zside][imax-1]+dz[zside],double(imax)*j_total/double(ncoils),32);
			}
		}
	}
	
}


