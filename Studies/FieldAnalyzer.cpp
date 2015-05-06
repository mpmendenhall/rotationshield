/* 
 * FieldAnalyzer.cpp, part of the RotationShield program
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

#include "FieldAnalyzer.hh"
#include "ProgressBar.hh"

double FieldAnalyzer::simpsonCoeff(unsigned int n, unsigned int ntot) {
    if(ntot == 1) return 1.0;
    if(ntot == 2) return 0.5;
    if(n==0 || n==ntot-1) return 1.0/(3.0*double(ntot-1));
    if(n%2==0) return 2.0/(3.0*double(ntot-1));
    return 4.0/(3.0*double(ntot-1));
}

void FieldAnalyzer::survey(vec3 ll, vec3 ur, unsigned int nX, unsigned int nY, unsigned int nZ, std::ostream& statsout, std::ostream& datsout) const {
    printf("Surveying resulting fields...\n");
    
    assert(nX>0 && nY>0 && nZ>0);
    vec3 xcenter = (ll+ur)*0.5;    // center of cell, around to produce polynomial fit
    
    unsigned int n[3] = { nX-1, nY-1, nZ-1 };
    vec3 dx = (ur - ll)/vec3(n[0],n[1],n[2]);
    for(int i=0; i<3; i++)
        if(n[i]==0) {
            dx[i] = 0;
            ll[i] = 0.5*(ll[i]+ur[i]);
    }
    
    
    
    // set up results tables
    gsl_matrix* coords = gsl_matrix_alloc(nX*nY*nZ,3);
    gsl_vector* bfield[3];
    for(int i=0; i<3; i++) bfield[i] = gsl_vector_alloc(nX*nY*nZ);
    
    ProgressBar pb(nX*nY);
    
    // scan over survey points
    vec3 bsXYZ(0,0,0);
    vec3 bssXYZ(0,0,0);
    double avgBx0 = 0;
    double avgBx1 = 0;
    int c=0; // incrementing coordinate index
    for(n[0] = 0; n[0] < nX; n[0]++) {
    
        vec3 bsYZ(0,0,0);
        vec3 bssYZ(0,0,0);
        for(n[1] = 0; n[1] < nY; n[1] ++) {
            
            pb.update(n[0]*nY+n[1]);
            
            vec3 bsZ(0,0,0);
            vec3 bssZ(0,0,0);
            for(n[2] = 0; n[2] < nZ; n[2]++)
            {
                vec3 x = ll + dx * vec3(n[0],n[1],n[2]);
                vec3 b = FS->fieldAt(x);
                for(int i=0; i<3; i++) gsl_vector_set(bfield[i],c,b[i]);
                
                bsZ += b*simpsonCoeff(n[2],nZ);
                bssZ += b*b*simpsonCoeff(n[2],nZ);
                
                for(int i=0; i<3; i++) datsout << std::setw(10) << std::setprecision(6) << x[i] << " ";
                datsout << "\t";
                for(int i=0; i<3; i++) datsout << std::setw(16) << std::setprecision(8) << gsl_vector_get(bfield[i],c) << " ";
                datsout << std::endl;
                
                for(int i=0; i<3; i++) gsl_matrix_set(coords, c, i, x[i] - xcenter[i]);
                
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
    
    
    // summary data and field fit
    double bRMS = sqrt(bssXYZ.sum());
    statsout << "#Field Shape" << std::endl;
    statsout << "#<B> = (" << bsXYZ[0] << " " << bsXYZ[1] << " " << bsXYZ[2] << ")\n";
    statsout << "#<B^2> = (" << bssXYZ[0] << " " << bssXYZ[1] << " " << bssXYZ[2] << ")\n";
    statsout << "#sigmaBx/Bx0 = " << sqrt(bssXYZ[0] - bsXYZ[0]*bsXYZ[0])/bsXYZ[0] << std::endl;
    statsout << "#<dBx/dx>/Bx0 = " << (avgBx1-avgBx0)/(bsXYZ[0]*(ur[0]-ll[0])) << std::endl;
    double resids;
    char axisnames[] = "xyz";
    
    unsigned int polOrder = 0;
    for(unsigned int i=1; i<=8; i++)
        if(nX>i && nY>i && nZ>i) polOrder = i;
    for(int i=0; i<3; i++) {
        Polynomial<3,double> p = Polynomial<3,double>::lowerTriangleTerms(polOrder);
        resids = polynomialFit(coords, bfield[i], p);
        p.prune(1e-9*bRMS);
        resids = polynomialFit(coords, bfield[i], p);
        gsl_vector_free(bfield[i]);
        statsout << "#" << polOrder << "th order polynomial (" << p.terms.size() << " terms) fit to B" << axisnames[i]
            << " centered at " << xcenter << ", rms residual = " << resids << std::endl;
        p.tableForm(statsout);
    }
    gsl_matrix_free(coords);
    
}


void FieldAnalyzer::visualizeSurvey(vec3 ll, vec3 ur, unsigned int nX, unsigned int nY, unsigned int nZ) const {
    
    if(nX*nY*nZ) {
        unsigned int n[3] = { nX-1, nY-1, nZ-1 };
        vec3 dx = (ur - ll)/vec3(n[0],n[1],n[2]);
        for(int i=0; i<3; i++) if(n[i]==0) { dx[i] = 0; ll[i] = 0.5*(ll[i]+ur[i]); }
        
        for(n[0] = 0; n[0] < nX; n[0]++) {
            for(n[1] = 0; n[1] < nY; n[1] ++) {
                for(n[2] = 0; n[2] < nZ; n[2]++) {
                    
                    vec3 x = ll+vec3(n[0],n[1],n[2])*dx;
                    vec3 b = FS->fieldAt(x);
                    std::cout << "\tx = " << x << "\t\tB = " << b << std::endl;
                    
                    float bmag = b.mag();
                    if(!bmag) continue;
                    b *= log(1+pow(bmag,0.25))/bmag;
                    
                    vsr::startRecording();
                    Line(x-b*0.02,x+b*0.02).visualizeDirected();
                    vsr::stopRecording();
                }
            }
        }
    } else {
        for(float i=0; i<nX+nY+nZ; i++) {
            double l = float(i)/(nX+nY+nZ-1);
            vec3 x = ll * (1-l) + ur * l;
            vec3 b = FS->fieldAt(x);
            std::cout << "\tx = " << x << "\t\tB = " << b << std::endl;
            
            float bmag = b.mag();
            if(!bmag) continue;
            b *= log(1+pow(bmag,0.25))/bmag;
            vsr::startRecording();
            Line(x-b*0.02,x+b*0.02).visualizeDirected();
            vsr::stopRecording();
        }
    }
}

