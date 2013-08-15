/// Short vertical coil
/// N=30, radius=1.m, length=2.m,
/// shield radius nominal 10cm outside coil; shield 40cm longer than coil
/// cell dimensions y=40cm, x=7.5cm, z=10cm, inner edge 5cm from center in x (field) direction
void shortCoil(std::ostream& gridf, std::ostream& fitf) {
	
	MixedSource* ms = new MixedSource();
	ms->retain();
	mdouble cradius = 1.;
	mdouble shieldDistance = 0.1;
	mdouble sradius = cradius + shieldDistance;
	mdouble clen = 2.;
	mdouble slen = clen + 0.40;
	
	// construct coil
	VarVec<mdouble> pert = VarVec<mdouble>(1);
	//pert[0] = -0.00295;
	pert[0] = 0;
	CosThetaBuilder b = CosThetaBuilder(15, cradius, clen, &shiftPositioner);
	b.regularCoil(*ms,&pert);
	ms->visualize();
	
	// construct shield
	FieldEstimator2D fe = FieldEstimator2D();
	fe.addsource(vec2(-clen/2.0,cradius),1.0);
	fe.addsource(vec2(clen/2.0,cradius),1.0);
	
	ShieldBuilder* G = new ShieldBuilder(128);
	G->makeOptCyl(10, 20, sradius, -slen/2, slen/2, new PlaneSource(Plane(),10000) );
	
	SymmetricSolver* sp = new SymmetricSolver(G);
	
	// solve shield for coil incident
	sp->solve();
	sp->calculateIncident(ms);
	sp->calculateResult();
	
	// map combined shield+coil data
	ms->addsource(sp);
	ms->visualize();
	
	FieldAnalyzer FA = FieldAnalyzer(ms);
	vec3 cell_ll(0.05,-0.05,-0.20);
	vec3 cell_ur(0.125,0.05,0.20);
	FA.visualizeSurvey(cell_ll,cell_ur,3,3,11);
	FA.survey(cell_ll,cell_ur,9,9,31,fitf,gridf);
	
	ms->release();
}