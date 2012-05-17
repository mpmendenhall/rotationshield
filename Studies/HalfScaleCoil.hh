#include <cassert>

/// Half-scale coil
/// N=30, radius=0.324, length=2.146,
/// distortion a = -0.00295
/// shield radius nominal 3.45cm outside coil
/// as-built radius delta 3.175cm?
void halfScaleCoil(mdouble shieldDistance, std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.324;
	mdouble sradius = cradius + shieldDistance;
	mdouble clen = 2.146;
	mdouble slen = clen + 0.20;
	
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
	
	//FA.survey(sampleLL*0.5,sampleUR*0.5,9,9,9,fitf,gridf);
	FA.survey(vec3(-.07,-.03,-.12),vec3(.07,.03,.12),15,5,5,fitf,gridf);
	
	ms->release();
}


/// Half-scale coil
/// N=30, radius=0.324, length=2.146,
/// distortion a = -0.00295
/// optimized shield radius nominal 3.45cm outside coil
/// as-built radius delta 3.175cm
void halfScaleTranslated(vec3 translation, std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.324;
	mdouble sradius = cradius + 0.03175;
	mdouble clen = 2.146;
	mdouble slen = clen + 0.20;
	
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;
	doubleTrans dt;
	dt.trans1 = translation;
	dt.trans2 = translation;
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner, &fancytrans);
	b.regularCoil(*ms,&pert,&dt);
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
	
	assert(false); // TODO
	//FA.survey(-sampleUR*0.5,sampleUR*0.5,11,11,11,fitf,gridf);
	
	ms->release();
}

void halfScaleTilted(vec3 translation, std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 0.324;
	mdouble sradius = cradius + 0.03175;
	mdouble clen = 2.146;
	mdouble slen = clen + 0.20;
	
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	pert[0] = -0.00295;
	doubleTrans dt;
	dt.trans1 = translation;
	dt.trans2 = -translation;
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner, &fancytrans);
	b.regularCoil(*ms,&pert,&dt);
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
	
	assert(false); // TODO
	//FA.survey(-sampleUR*0.5,sampleUR*0.5,11,11,11,fitf,gridf);
	
	
	ms->release();
}