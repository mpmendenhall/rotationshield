#include "Studies.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include <cassert>

//--

coilShield::coilShield(): length(0), radius(0), mu(10000.), cSegs(10), vSegs(20), pSegs(128) {
	for(unsigned int i=0; i<2; i++) {
		endcap_delta_z[i] = endcap_delta_or[i] = endcap_ir[i] = 0;
		endcap_mu[i] = 10000;
		eSegs[i] = 0;
	}
}

Stringmap coilShield::getInfo() const {
	Stringmap m;
	m.insert("length",length);
	m.insert("radius",radius);
	m.insert("cSegs",itos(cSegs));
	m.insert("vSegs",itos(vSegs));
	m.insert("pSegs",itos(pSegs));
	m.insert("mu",mu);
	for(int i=0; i<2; i++) {
		if(eSegs[i]) {
			m.insert("endcap_segs_"+itos(i),itos(eSegs[i]));
			m.insert("endcap_ir_"+itos(i),endcap_ir[i]);
			m.insert("endcap_delta_or_"+itos(i),endcap_delta_or[i]);
			m.insert("endcap_delta_z_"+itos(i),endcap_delta_z[i]);
			m.insert("endcap_mu_"+itos(i),endcap_mu[i]);
		}
	}
	return m;
}

void coilShield::construct(MixedSource& ms, CosThetaBuilder* ct) const {
	
	FieldEstimator2Dfrom3D fe(&ms);
	
	ShieldBuilder* G = new ShieldBuilder(pSegs);
	G->makeOptCyl(cSegs, vSegs, radius, -length/2., length/2., new PlaneSource(Plane(),mu), &fe);
	
	for(int i=0; i<2; i++) {
		if(!eSegs[i]) continue;
		int z = 2*i-1;
		float halfz = z*(length + endcap_delta_z[i])/2;
		G->OptCone(eSegs[i]/3, eSegs[i]-(eSegs[i]/3),
					vec2(halfz, z<0? endcap_ir[i] : radius+endcap_delta_or[i]),
					vec2(halfz, z<0? radius+endcap_delta_or[i] : endcap_ir[i]),
					new PlaneSource(Plane(),endcap_mu[i]), &fe);
	}
	
	SymmetricSolver* sp = new SymmetricSolver(G);
	
	sp->solve();
	sp->calculateIncident(&ms);
	sp->calculateResult();
	
	ms.addsource(sp);
}

//--

Stringmap fieldCell::getInfo() const {
	Stringmap m;
	m.insert(std::string("ll"),vtos(vec2doublevec<3,mdouble>(ll)));
	m.insert(std::string("ur"),vtos(vec2doublevec<3,mdouble>(ur)));
	m.insert("nx",nx);
	m.insert("ny",ny);
	m.insert("nz",nz);
	return m;
}

//--

nEDM_Geom::nEDM_Geom(const std::string& dir): coil(NULL), shield(NULL), cell(NULL),
saveGrid(true), basedir(dir), ms(NULL) { }

nEDM_Geom::~nEDM_Geom() {
	if(ms) ms->release();
}
		
void nEDM_Geom::construct() {
	if(ms) ms->release();
	ms = new MixedSource();
	ms->retain();
	assert(coil);
	coil->myCap = CosThetaBuilder::CAP_ARC;
	coil->buildCoil(*ms);
	ms->visualize();
	if(shield) {
		shield->construct(*ms,coil);
		ms->visualize();
	}
}

void nEDM_Geom::takeSample(const std::string& sName) {
	assert(cell);
	if(!ms) {
		printf("No field sources specified!\n");
		return;
	}
	
	// set up files
	std::string fieldspath = basedir+"/"+sName;
	makePath(fieldspath);
	std::ofstream fieldsout;
	std::ofstream statsout;
	std::ostream nullout(NULL);
	if(saveGrid) fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
	statsout.open((fieldspath+"/Fieldstats.txt").c_str());
	
	// save geometry summary info
	QFile qOut;
	if(coil)
		coil->writeInfo(qOut);
	if(shield)
		qOut.insert("shield",shield->getInfo());
	if(cell)
		qOut.insert("cell",cell->getInfo());
	qOut.commit(fieldspath+"/GeomInfo.txt");
	
	// run analyzer
	FieldAnalyzer FA = FieldAnalyzer(ms);
	FA.visualizeSurvey(cell->ll,cell->ur,cell->vx,cell->vy,cell->vz);
	FA.survey(cell->ll,cell->ur,cell->nx,cell->ny,cell->nz,statsout,saveGrid?fieldsout:nullout);

	// cleanup
	if(saveGrid) fieldsout.close();
	statsout.close();
}





//----------------------------------------------------------------
// Old studies... to be re-written with newer utility functions
//----------------------------------------------------------------

/*

/// Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
/// shield radius nominal 6.9cm outside coil
/// shield 40cm longer than coil
void fullScaleCoil(std::ostream& gridf, std::ostream& fitf);

/// Bare Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
void bareFullScaleCoil(std::ostream& gridf, std::ostream& fitf);

void fullScaleCoil(std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.648; //64.8cm radius
	mdouble shieldDistance = 0.069; //6.9cm spacing to shield
	mdouble sradius = cradius + shieldDistance;
	mdouble clen = 4.292;
	mdouble slen = clen + 0.40;
	
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner);
	b.regularCoil(*ms,&pert);
	ms->visualize();
	
	FieldEstimator2D fe = FieldEstimator2D();
	fe.addsource(vec2(-clen/2.0,cradius),1.0);
	fe.addsource(vec2(clen/2.0,cradius),1.0);
	
	ShieldBuilder* G = new ShieldBuilder(128);
	G->makeOptCyl(10, 20, sradius, -slen/2, slen/2, new PlaneSource(Plane(),10000) );
	
	SymmetricSolver* sp = new SymmetricSolver(G);
	
	sp->solve();
	sp->calculateIncident(ms);
	sp->calculateResult();
	
	ms->addsource(sp);
	ms->visualize();
	
	FieldAnalyzer FA = FieldAnalyzer(ms);
	FA.visualizeSurvey(vec3(-0.25,-0.25,0.0),vec3(0.25,0.25,3.0),3,3,11);
	FA.survey(vec3(-0.25,-0.25,0.0),vec3(0.25,0.25,3.0),9,9,31,fitf,gridf);
	
	ms->release();
}

void bareFullScaleCoil(std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.648;
	mdouble clen = 4.292;
	
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner);
	b.regularCoil(*ms,&pert);
	ms->visualize();
	
	FieldEstimator2D fe = FieldEstimator2D();
	fe.addsource(vec2(-clen/2.0,cradius),1.0);
	fe.addsource(vec2(clen/2.0,cradius),1.0);
	
	FieldAnalyzer FA = FieldAnalyzer(ms);
	FA.visualizeSurvey(vec3(-0.25,-0.25,0.0),vec3(0.25,0.25,3.0),3,3,11);
	//FA.survey(vec3(-0.40,-0.40,-3.0),vec3(0.40,0.40,3.0),81,81,601,fitf,gridf);
	FA.survey(vec3(-0.40,-0.40,-3.0),vec3(0.40,0.40,3.0),21,21,151,fitf,gridf);
	
	ms->release();
}


/// full example of doing everything the hard way...
void mi_sampleShieldVis(std::deque<std::string>&, std::stack<std::string>&) {
	// set up output paths
	std::string basepath = "../";
	std::string projectpath = basepath+"/CosThetaCoil/";
	makePath(projectpath);
	std::string fieldspath = projectpath+"/Fields/";
	makePath(fieldspath);
	std::string gfpath = projectpath+"/GFs/";
	makePath(gfpath);
	
	// set up output files for field scans
	std::ofstream fieldsout;
	std::ofstream statsout;
	fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
	statsout.open((fieldspath+"/Fieldstats.txt").c_str());
	
	// set up cos theta coil
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.648;	// coil radius
	mdouble clen = 4.292;		// coil length
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;	// coil distortion parameter
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner);
	b.regularCoil(*ms,&pert);
	ms->visualize();
	//vsr::pause();
	
	// optional: construct shield
	if(1) {
		mdouble shieldDistance = 0.069;				// distance of shield from wires, 6.9cm
		mdouble sradius = cradius + shieldDistance;	// shield radius
		mdouble slen = clen + 0.40;					// shield length, coil + 40cm
		FieldEstimator2D* fe = new FieldEstimator2D();
		fe->addsource(vec2(-clen/2.0,cradius),1.0);
		fe->addsource(vec2(clen/2.0,cradius),1.0);
		ShieldBuilder* G = new ShieldBuilder(128);	// 128 segments around shield builder
		G->makeOptCyl(10, 20, sradius, -slen/2, slen/2, new PlaneSource(Plane(),10000), fe); // optimized grid with 10+20 divisions along z
		// solve for shield response function
		SymmetricSolver* sp = new SymmetricSolver(G);
		sp->solve();
		// calculate response to coil fields
		sp->calculateIncident(ms);
		sp->calculateResult();
		ms->addsource(sp);
		ms->visualize();
	}
	
	// example of probing field value
	vec3 b0 = ms->fieldAt(vec3(0,0,0));
	printf("Field at center: ");
	b0.display();
	
	// survey data points on grid
	FieldAnalyzer FA = FieldAnalyzer(ms);
	//        lower left corner,     upper right corner, nx,ny,nz, output files
	FA.survey(vec3(0.0,0.0,0.0),vec3(0.20,0.20,0.5),10,10,10,statsout,fieldsout);
	
	
	// cleanup
	ms->release();
	fieldsout.close();
	statsout.close();
	vsr::pause();
}

*/
