/// Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
/// shield radius nominal 6.9cm outside coil
/// shield 40cm longer than coil
void fullScaleCoil(std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.648;
	mdouble shieldDistance = 0.069;
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

/// Bare Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
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