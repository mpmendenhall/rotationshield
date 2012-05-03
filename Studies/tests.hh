/// \file "tests.hh" \brief Tests to make sure everything is working properly
#ifndef TESTS_HH
/// Make sure this header is only loaded once
#define TESTS_HH 1

#include "SymmetricSolver.hh"
#include "GenericSolver.hh"
#include "FieldSource.hh"
#include "Boxel.hh"
#include "ShieldBuilder.hh"
#include "analysis.hh"
#include "CosThetaBuilder.hh"

/// Print a comparison between expected and calculated values for a data point
bool compareResults(mdouble a, mdouble b, const char* label = 0x0) {
	bool pass = true;
	if(label) printf("%s:\n",label);
	printf("\tTest Result:\t%.14e\n\t  Reference:\t%.14e\n\t\t   error:\t%.1f%%",double(a),double(b),double(100*(a-b)/a));
	if(fabs(100*(a-b)/a) > 5.0 ) { printf(" <<-------- ******* TEST FAILED ******"); pass = false; }
	printf("\n");
	fflush(stdout);
	return pass;
}

/// Basic test that everything is working well
bool reference_sanity_check()
{
	MixedSource* MxS = new MixedSource();
	MxS->retain();
	
	CosThetaBuilder(17, 0.61, 3.92).regularCoil(*MxS);
	MxS->visualize();
	
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan = vec3(0.15,0,0);
	vec3 yscan = vec3(0,0.06,0);
	vec3 zscan = vec3(0,0,0.25);
	
	bool pass = true;
	mdouble b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,8.53210605848347e-01,"coil origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,8.53217812021647e-01-8.53210605848347e-01,"coil +x edge");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,8.53163006651345e-01-8.53210605848347e-01,"coil +y edge");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,8.54088402519457e-01-8.53210605848347e-01,"coil +z edge");
	pass &= compareResults(MxS->fieldAt(xscan*-1.0)[0]-b0,8.53217812021647e-01-8.53210605848347e-01,"coil -x edge");
	pass &= compareResults(MxS->fieldAt(yscan*-1.0)[0]-b0,8.53163006651345e-01-8.53210605848347e-01,"coil -y edge");
	pass &= compareResults(MxS->fieldAt(zscan*-1.0)[0]-b0,8.54088402519457e-01-8.53210605848347e-01,"coil -z edge");
	
	if(!pass) { vsr::Visr::W->pause(); return pass; }
	
	
	FieldEstimator2D fe;
	fe.addsource(vec2(-3.92/2,0.61),1.0);
	fe.addsource(vec2(3.92/2,0.61),1.0);
	
	ShieldBuilder SB = ShieldBuilder(128);
	//SB.makeOptCyl(20, 30, .6223, -3.9624/2, 3.9624/2, new Boxel(annulusSpec(),1e5,0.001), &fe);
	SB.makeOptCyl(20, 30, .6223, -3.9624/2, 3.9624/2, new PlaneSource(Plane(),10000.0), &fe);

	SymmetricSolver* s = new SymmetricSolver(&SB);
	s->retain();
	s->visualize();
	
	s->calculateIncident(MxS);
	s->visualize();
	
	s->solve();

	s->calculateResult();
	
	MxS->addsource(s);
	b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.60612752621056,"all origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,4.93544984829120e-04,"all +x edge");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-6.03581955409549e-05,"all +y edge");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-2.06183474346311e-04,"all +z edge");
	pass &= compareResults(MxS->fieldAt(xscan*-1.0)[0]-b0,4.93544984829342e-04,"all -x edge");
	pass &= compareResults(MxS->fieldAt(yscan*-1.0)[0]-b0,-6.03581955422872e-05,"all -y edge");
	pass &= compareResults(MxS->fieldAt(zscan*-1.0)[0]-b0,-2.06183474332544e-04,"all -z edge");
	
	vsr::Visr::W->pause();
	MxS->visualize();
	vsr::Visr::W->pause();
	
	s->release();
	MxS->release();
	
	return pass;
}


/// Faster test with a very crudely gridded shield
bool reference_simpleshield() {
	
	MixedSource* MxS = new MixedSource();
	CosThetaBuilder b = CosThetaBuilder(5, 0.55, 3.92);
	b.regularCoil(*MxS);
	FieldAnalyzer FA = FieldAnalyzer(MxS);
	
	FieldEstimator2D fe;
	fe.addsource(vec2(-3.92/2,0.55),1.0);
	fe.addsource(vec2(3.92/2,0.55),1.0);
	
	ShieldBuilder SB = ShieldBuilder(32);
	SB.makeOptCyl(10, 2, .6223, -3.9624/2, 3.9624/2, new PlaneSource(Plane(),10000.0), &fe);
	//SB.makeOptCyl(10, 2, .6223, -3.9624/2, 3.9624/2, new Boxel(annulusSpec(),100000.0,0.001), &fe);

	SymmetricSolver* s = new SymmetricSolver(&SB);
	s->retain();
	
	s->calculateIncident(MxS);
	s->visualize();
	vsr::Visr::W->pause();
	
	s->solve();
	s->calculateResult();
	
	MxS->addsource(s);
	MxS->visualize();
	
	//center scan lines
	vec3 origin(0,0,0);
	vec3 xscan(0.15,0,0);
	vec3 yscan(0,0.06,0);
	vec3 zscan(0,0,0.25);
	
	bool pass = true;
	mdouble b0 = MxS->fieldAt(origin)[0];
	pass &= compareResults(b0,1.61484131308468e+00,"origin");
	pass &= compareResults(MxS->fieldAt(xscan)[0]-b0,2.64365094181551e-03,"+x");
	pass &= compareResults(MxS->fieldAt(yscan)[0]-b0,-3.81123880558443e-04,"+y");
	pass &= compareResults(MxS->fieldAt(zscan)[0]-b0,-3.01945405914772e-04,"+z");
	pass &= compareResults(MxS->fieldAt(-xscan)[0]-b0,2.64365094181596e-03,"-x");
	pass &= compareResults(MxS->fieldAt(-yscan)[0]-b0,-3.81123880558221e-04,"-y");
	pass &= compareResults(MxS->fieldAt(-zscan)[0]-b0,-3.01945405912996e-04,"-z");
	
	vsr::Visr::W->pause();
	s->release();
	return pass;
}




#endif
