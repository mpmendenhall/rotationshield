#include "analysis.hh"


mdouble FieldAnalyzer::simpsonCoeff(int n, int ntot) {
	if(ntot == 1) return 1.0;
	if(ntot == 2) return 0.5;
	if(n==0 || n==ntot-1) return 1.0/(3.0*mdouble(ntot-1));
	if(n%2==0) return 2.0/(3.0*mdouble(ntot-1));
	return 4.0/(3.0*mdouble(ntot-1));
}

void FieldAnalyzer::survey(vec3 ll, vec3 ur, int nX, int nY, int nZ, std::ostream& statsout, std::ostream& datsout) const {	
	printf("Analyzing resultant fields... ");
	vec3 n = vec3(nX-1,nY-1,nZ-1);
	vec3 dx = (ur - ll)/n;
	for(int i=0; i<3; i++) if(n[i]==0) { dx[i] = 0; ll[i] = 0.5*(ll[i]+ur[i]); }
	vec3 b, x, bsZ, bsYZ, bsXYZ, bssZ, bssYZ, bssXYZ;
	mdouble avgBx0 = 0;
	mdouble avgBx1 = 0;
	
	// set up results tables
	gsl_matrix* coords = gsl_matrix_alloc(nX*nY*nZ,3);
	gsl_vector* bfield[3];
	for(int i=0; i<3; i++) bfield[i] = gsl_vector_alloc(nX*nY*nZ);
	
	// scan over survey points
	bsXYZ = vec3(); bssXYZ = vec3();
	int c=0;
	for(n[0] = 0; n[0] < nX; n[0]++) {
		bsYZ = vec3(); bssYZ = vec3();
		for(n[1] = 0; n[1] < nY; n[1] ++) {
			
			printf("*"); fflush(stdout);
			
			bsZ = vec3(); bssZ = vec3();
			for(n[2] = 0; n[2] < nZ; n[2]++)
			{
				x = ll+n*dx;
				b = FS->fieldAt(x);
				for(int i=0; i<3; i++) gsl_vector_set(bfield[i],c,b[i]);
				
				bsZ += b*simpsonCoeff(n[2],nZ);
				bssZ += b*b*simpsonCoeff(n[2],nZ);
				
				for(int i=0; i<3; i++) datsout << std::setw(16) << x[i] << " ";
				for(int i=0; i<3; i++) datsout << std::setw(20) << std::setprecision(12) << gsl_vector_get(bfield[i],c) << " ";
				datsout << std::endl;
				
				for(int i=0; i<3; i++) gsl_matrix_set(coords,c,i,x[i]);
				
				c++;
			}
			bsYZ += bsZ*simpsonCoeff(n[1],nY);
			bssYZ += bssZ*simpsonCoeff(n[1],nY);
		}
		if(n[0]==0) avgBx0 = bsYZ[0];
		if(n[0]==nX-1) avgBx1 = bsYZ[0];
		bsXYZ += bsYZ*simpsonCoeff(n[0],nX);
		bssXYZ += bssYZ*simpsonCoeff(n[0],nX);
	}
	printf("\n");
	
	
	// summary data and field fit
	mdouble bRMS = sqrt(bssXYZ.sum());
	statsout << "#Field Shape" << std::endl;
	statsout << "#<B> = (" << bsXYZ[0] << " " << bsXYZ[1] << " " << bsXYZ[2] << ")\n";
	statsout << "#<B^2> = (" << bssXYZ[0] << " " << bssXYZ[1] << " " << bssXYZ[2] << ")\n";
	statsout << "#sigmaBx/Bx0 = " << sqrt(bssXYZ[0] - bsXYZ[0]*bsXYZ[0])/bsXYZ[0] << std::endl;
	statsout << "#<dBx/dx>/Bx0 = " << (avgBx1-avgBx0)/(bsXYZ[0]*(ur[0]-ll[0])) << std::endl;
	mdouble resids;
	char axisnames[] = "xyz";
	
	unsigned int polOrder = 4;
	if(nX>5 && nY>5 && nZ>5) polOrder = 5;
	if(nX>6 && nY>6 && nZ>6) polOrder = 6;
	Polynomial<3,mdouble> p = Polynomial<3,mdouble>::lowerTriangleTerms(polOrder);
	for(int i=0; i<3; i++) {
		resids = polynomialFit(coords, bfield[i], p);
		gsl_vector_free(bfield[i]);
		statsout << "#" << polOrder << "th order polynomial (" << p.terms.size() << " terms) fit to B" << axisnames[i] << " centered at < 0 0 0 >, rms residual = " << resids << std::endl;
		Polynomial<3,mdouble> p2 = p;
		p2.prune(1e-9*bRMS);
		p2.tableForm(statsout);
		//statsout << "#4th order Polynomial for B" << axisnames[i] << " recentered at " << (ll+ur)*0.5 << std::endl;
		//p2 = p.recentered((ll+ur)*0.5);
		//p2.prune(1e-9);
		//p2.tableForm(statsout);
	}
	gsl_matrix_free(coords);
	
}


void FieldAnalyzer::visualizeSurvey(vec3 ll, vec3 ur, int nX, int nY, int nZ) const {	
	
	vec3 n = vec3(nX-1,nY-1,nZ-1);
	vec3 dx = (ur - ll)/n;
	vec3 x,b;
	for(int i=0; i<3; i++) if(n[i]==0) { dx[i] = 0; ll[i] = 0.5*(ll[i]+ur[i]); }
	
	vsr::startRecording(); 
	vsr::clearWindow();
	FS->visualize(false);
	
	for(n[0] = 0; n[0] < nX; n[0]++) {
		for(n[1] = 0; n[1] < nY; n[1] ++) {
			for(n[2] = 0; n[2] < nZ; n[2]++) {
				
				x = ll+n*dx;
				b = FS->fieldAt(x);
				float bmag = b.mag();
				if(!bmag) continue;
				b *= log(1+pow(bmag,0.25))/bmag;
				Line(x-b*0.02,x+b*0.02).visualizeDirected();
			}
		}
	}
	
	vsr::stopRecording();
}

