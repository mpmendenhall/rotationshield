/* 
 * Studies.cpp, part of the RotationShield program
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

#include "Studies.hh"
#include "strutils.hh"
#include "PathUtils.hh"
#include "FieldAdaptiveSurface.hh"
#include "LineSource.hh"
#include "UniformField.hh"
#include "SurfaceGeometry.hh"
#include "SurfaceProfiles.hh"
#include "SurfaceCurrentRS.hh"
#include "HoleDipolePerturbation.hh"
#include "GenericSolver.hh"
#include "MultiQuilibrator.hh"
#include "SymmetrizedSource.hh"
#include <cassert>

//--

Stringmap fieldCell::getInfo() const {
    Stringmap m;
    m.insert(string("ll"),vtos(vec2doublevec<3,double>(ll)));
    m.insert(string("ur"),vtos(vec2doublevec<3,double>(ur)));
    m.insert("nx",nx);
    m.insert("ny",ny);
    m.insert("nz",nz);
    return m;
}

//--

//-------------------------------------------
// Menu-driven user interface actions
//-------------------------------------------

void mi_setFCrange(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float urz = S->popFloat();
    float ury = S->popFloat();
    float urx = S->popFloat();
    float llz = S->popFloat();
    float lly = S->popFloat();
    float llx = S->popFloat();
    SC->cell.ll = vec3(llx,lly,llz);
    SC->cell.ur = vec3(urx,ury,urz);
}

void mi_setFCgrid(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->cell.nz = S->popInt();
    SC->cell.ny = S->popInt();
    SC->cell.nx = S->popInt();
}

void mi_setSaveGrid(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->cell.saveGrid = S->popInt();
}

void mi_outDir(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->basedir += "/" + S->popString();
}

//-------------------------------------------

void mi_toggleUpdateVis(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->updateVis = !SC->updateVis;
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addLineCurrent(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float j = S->popFloat();
    float ez = S->popFloat();
    float ey = S->popFloat();
    float ex = S->popFloat();
    float sz = S->popFloat();
    float sy = S->popFloat();
    float sx = S->popFloat();
    
    SC->IncidentSource->addsource(new LineSource(vec3(sx,sy,sz), vec3(ex,ey,ez), j));
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addRingCurrent(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float j = S->popFloat();
    float r = S->popFloat();
    float z = S->popFloat();
    
    SC->IncidentSource->loop(z,r,j);
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_symmetrizeField(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    int p = S->popInt();
    if(!p) return;
    
    MixedSource* MS = new MixedSource();
    MS->addsources(SC->IncidentSource);
    SymmetrizedSource* SS = new SymmetrizedSource(MS,p>0);
    SC->IncidentSource->clear();
    SC->IncidentSource->addsource(SS);
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addUnifB(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float bz = S->popFloat();
    float by = S->popFloat();
    float bx = S->popFloat();
    SC->IncidentSource->addsource(new UniformField(vec3(bx,by,bz)));
}

void mi_ClearIncident(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->IncidentSource->clear();
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_buildCosThetaExit(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->CTB.buildCoil(*SC->IncidentSource);
    S->mydeque->push_front(NameSelector::exit_control);
    if(SC->updateVis) SC->TotalField->visualize();
}

//---------------------------------------------

void mi_Solve(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    string savenm = S->popString();
    if(savenm=="none") savenm = "";
    SC->solve(savenm);
    SC->calculate_result();
    SC->TotalField->visualize();
}

void mi_Recalc(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    if(!SC->SS || !SC->RSC) { printf("Solver not yet generated; nothing done.\n"); return; }
    SC->calculate_result();
    SC->TotalField->visualize();
}

void mi_zeroResponse(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    if(!SC->SS || !SC->RSC) { printf("Solver not yet generated; nothing done.\n"); return; }
    SC->RSC->setFinalState(mvec(SC->RSC->nDF()));
    SC->TotalField->visualize();
}

void mi_meas(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    if(!SC->RSC) printf("Measuring noninteracting applied fields:\n");
    else printf("Measuring resulting fields:\n");
    SC->measureFields();
    SC->writeInfo();
    printf("Data collection complete.\n");
}

void mi_qSurvey(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    unsigned int i = S->popInt();
    if(Visualizable::vis_on) {
        SC->TotalField->visualize();
        SC->FA.visualizeSurvey(SC->cell.ll, SC->cell.ur, i, 0, 0);
    }
}

void mi_addSingular(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float c = S->popFloat();
    unsigned int i = S->popInt();
    if(!SC->SS || !SC->RSC) { printf("Solver not yet generated; nothing done.\n"); return; }
    if(i >= SC->RSC->nDF()) { printf("Maximum vector number is %i; nothing done.\n", SC->RSC->nDF()-1); return; }
    SC->add_singular_state(i,c);
    SC->TotalField->visualize();
}

void mi_setSingularEpsilon(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float ep = S->popFloat();
    if(!SC->SS) { printf("Solver not yet generated; nothing done.\n"); return; }
    SC->SS->set_singular_epsilon(ep);
    mi_Recalc(S);
}

//---------------------------------------------

void mi_setPhi(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    SC->initReactiveSet(S->popInt());
}

void mi_addTube(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float o = S->popFloat();
    unsigned int ns = S->popInt();
    float mu = S->popFloat();
    float r = S->popFloat();
    float er = S->popFloat();
    float ez = S->popFloat();
    float sr = S->popFloat();
    float sz = S->popFloat();
    
    RoundedTube* RT = new RoundedTube(vec2(sz,sr), vec2(ez,er), r);
    DVFunc1<2,double>* FAS = SC->adaptSurface(RT,o);
    CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
    SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, SC->RSC->nPhi, ns);
    RS->setSurfaceResponse(SurfaceI_Response(mu));
    RS->point_ixn = false;
    SC->RSC->addSet(RS);
    
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addSlab(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float o = S->popFloat();
    unsigned int ns = S->popInt();
    float mu = S->popFloat();
    float er = S->popFloat();
    float r = S->popFloat();
    float z = S->popFloat();
    
    RoundedSlab* SB = new RoundedSlab(z,r,2*er);
    DVFunc1<2,double>* FAS = SC->adaptSurface(SB,o);
    CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
    SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, SC->RSC->nPhi, ns);
    RS->setSurfaceResponse(SurfaceI_Response(mu));
    RS->point_ixn = false;
    SC->RSC->addSet(RS);
    
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addWiggleSlab(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float o = S->popFloat();
    unsigned int ns = S->popInt();
    float mu = S->popFloat();
    float a = S->popFloat();
    unsigned int nw = (unsigned int)S->popInt();
    float er = S->popFloat();
    float r = S->popFloat();
    float z = S->popFloat();
    
    WiggleSlab* SB = new WiggleSlab(z, r, 2*er, nw, a);
    DVFunc1<2,double>* FAS = SC->adaptSurface(SB,o);
    CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
    SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, SC->RSC->nPhi, ns);
    RS->setSurfaceResponse(SurfaceI_Response(mu));
    RS->point_ixn = false;
    SC->RSC->addSet(RS);
    
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addBall(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float o = S->popFloat();
    unsigned int ns = S->popInt();
    float mu = S->popFloat();
    float r = S->popFloat();
    float z = S->popFloat();
    
    Arc2D* B = new Arc2D(r, -M_PI, 0, vec2(z,0));
    DVFunc1<2,double>* FAS = SC->adaptSurface(B,o);
    CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
    SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, SC->RSC->nPhi, ns);
    RS->setSurfaceResponse(SurfaceI_Response(mu));
    RS->point_ixn = false;
    SC->RSC->addSet(RS);
    
    if(SC->updateVis) SC->TotalField->visualize();
}

void mi_addTorus(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float o = S->popFloat();
    unsigned int ns = S->popInt();
    float mu = S->popFloat();
    float r2 = S->popFloat();
    float r1 = S->popFloat();
    float z = S->popFloat();
    if(fabs(r2)>fabs(r1)) std::swap(r1,r2);
    
    Arc2D* B = new Arc2D(r2, -M_PI, M_PI, vec2(z,r1));
    DVFunc1<2,double>* FAS = SC->adaptSurface(B,o);
    CylSurfaceGeometry* SG = new CylSurfaceGeometry(FAS);
    SurfaceCurrentRS* RS = new SurfaceCurrentRS(SG, SC->RSC->nPhi, ns);
    RS->setSurfaceResponse(SurfaceI_Response(mu));
    RS->point_ixn = false;
    SC->RSC->addSet(RS);
    
    if(SC->updateVis) SC->TotalField->visualize();
}

//---------------------------------------------

void mi_punchHole(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);
    float r = S->popFloat();
    float z = S->popFloat();
    float y = S->popFloat();
    float x = S->popFloat();
    
    const vector<ReactiveSet*>& v = SC->RSC->getSets();
    if(!v.size()) {
        std::cout << "No previous surface defined to puncture!\n";
        return;
    }
    SurfaceCurrentRS* RS = dynamic_cast<SurfaceCurrentRS*>(v.back());
    
    double d2;
    vec2 l = RS->mySurface->closestPoint(vec3(x,y,z),d2);
    std::cout << "Adding hole at " << (*RS->mySurface)(l) << " of radius " << r << "\n";
    
    HoleDipolePerturbation* h = new HoleDipolePerturbation(*RS->mySurface, l, r);
    h->hide_ixn.insert(RS->RS_UID);
    SC->PTB->addSet(h);
}

void mi_equilibratePtb(StreamInteractor* S) {
    SystemConfiguration* SC = dynamic_cast<SystemConfiguration*>(S);

    // solve perturbation system
    SC->PTB->incidentState = SC->PTB->getFullReactionTo(SC->IncidentSource);
    GenericSolver GS;
    GS.solve(*SC->PTB);
    
    // equilibrate system
    MultiQuilibrator MQ;
    MQ.addSet(SC->RSC, SC->SS);
    MQ.addSet(SC->PTB, &GS);
    while(MQ.step() > 1e-6)
        if(SC->updateVis)
            SC->TotalField->visualize();
}

//---------------------------------------------




//-------------------------------------------
//
//-------------------------------------------

SystemConfiguration::SystemConfiguration():
basedir(getEnvSafe("ROTSHIELD_OUT","./")),
RSC(new MagRSCombiner(32)),
PTB(new MagRSCombiner(1)),
IncidentSource(new MagExtField()),
TotalField(new MixedSource()),
FA(TotalField),
SS(NULL),
outDir("Output directory", &mi_outDir, this),
setFCrange("Set Measurement Range", &mi_setFCrange, this),
setFCgrid("Set Measurement Gridding", &mi_setFCgrid, this),
setSaveGrid("Enable/disable gridded output", &mi_setSaveGrid, this),
OMcell("Measurement Cell Options"),
addLineCurrent("Add line segment current", &mi_addLineCurrent, this),
addUnifB("Add uniform field", &mi_addUnifB, this),
clearIncident("Remove all field sources", &mi_ClearIncident, this),
buildCosThetaExit("Build coil and exit", &mi_buildCosThetaExit, this),
addRingCurrent("Add ring current", &mi_addRingCurrent, this),
symmetrizeField("(Anti)Symmetrize field", &mi_symmetrizeField, this),
toggleUpdateVis("Toggle visualization updating", &mi_toggleUpdateVis, this),
OMfieldsrc("Field Source Options"),
setPhi("Set rotational symmetry grid", &mi_setPhi, this),
addSlab("Add circular slab", &mi_addSlab, this),
addTube("Add tube", &mi_addTube, this),
addBall("Add ball", &mi_addBall, this),
addTorus("Add torus", &mi_addTorus, this),
addWiggleSlab("Add wiggly slab", &mi_addWiggleSlab, this),
OMsurfaces("Boundary conditions"),
doSolve("Solve boundary condition interactions",&mi_Solve,this),
doApply("(Re)apply incident field to boundaries",&mi_Recalc,this),
zeroResponse("Zero out surface response state",&mi_zeroResponse,this),
doMeas("Take field measurement",&mi_meas,this),
qSurvey("Check field along cell diagonal",&mi_qSurvey,this),
addSingular("Add enumerated singular state", &mi_addSingular, this),
setSingularEpsilon("Set SVD pseudo-inverse threshold", &mi_setSingularEpsilon, this),
punchHole("Add hole perturbation to superconducting surface", &mi_punchHole, this),
equilibratePtb("Equilibrate field perturbations", &mi_equilibratePtb, this) {
    
    IncidentSource->retain();
    TotalField->retain();
    RSC->retain();
    
    TotalField->addsource(IncidentSource);
    TotalField->addsource(RSC);
    TotalField->addsource(PTB);
    
    outDir.addArg("rel. to $ROTSHIELD_OUT",".");
    //
    setFCrange.addArg("x min","-0.20");
    setFCrange.addArg("y min","-0.20");
    setFCrange.addArg("z min","-0.20");
    setFCrange.addArg("x max","0.20");
    setFCrange.addArg("y max","0.20");
    setFCrange.addArg("z max","0.20");
    //
    setFCgrid.addArg("nx","5");
    setFCgrid.addArg("ny","5");
    setFCgrid.addArg("nz","5");
    //
    setSaveGrid.addArg("Enable","1");
    //
    OMcell.addChoice(&setFCrange,"range");
    OMcell.addChoice(&setFCgrid,"grid");
    OMcell.addChoice(&setSaveGrid,"svgrd");
    OMcell.addChoice(&InputRequester::exitMenu,"x");
    
    addLineCurrent.addArg("x0","0");
    addLineCurrent.addArg("y0","0");
    addLineCurrent.addArg("z0","-1");
    addLineCurrent.addArg("x1","0");
    addLineCurrent.addArg("y1","0");
    addLineCurrent.addArg("z1","1");
    addLineCurrent.addArg("I","1");
    //
    addUnifB.addArg("Bx","1");
    addUnifB.addArg("By","1");
    addUnifB.addArg("Bz","1");
    //
    addRingCurrent.addArg("z","0");
    addRingCurrent.addArg("r","1");
    addRingCurrent.addArg("I","1");
    //
    symmetrizeField.addArg("p","1");
    //
    OMfieldsrc.addChoice(&addLineCurrent,"l");
    OMfieldsrc.addChoice(&addRingCurrent,"r");
    OMfieldsrc.addChoice(&addUnifB,"u");
    CTB.OMcoil.addChoice(&buildCosThetaExit,"x");
    OMfieldsrc.addChoice(&CTB.OMcoil,"c");
    OMfieldsrc.addChoice(&clearIncident,"d");
    OMfieldsrc.addChoice(&symmetrizeField,"s");
    OMfieldsrc.addChoice(&InputRequester::exitMenu,"x");
    
    setPhi.addArg("nPhi", std::to_string(RSC->nPhi));
    //
    addSlab.addArg("z","-0.7");
    addSlab.addArg("r","0.5");
    addSlab.addArg("end radius","0.05");
    addSlab.addArg("mu","0");
    addSlab.addArg("nZ","11");
    addSlab.addArg("adapt","0");
    //
    addTube.addArg("z0","-0.6");
    addTube.addArg("r0","0.6");
    addTube.addArg("z1","0.6");
    addTube.addArg("r1","0.6");
    addTube.addArg("end radius","0.05");
    addTube.addArg("mu","10000");
    addTube.addArg("nZ","17");
    addTube.addArg("adapt","0");
    //
    addBall.addArg("z","0");
    addBall.addArg("r","0.5");
    addBall.addArg("mu","0");
    addBall.addArg("nZ","11");
    addBall.addArg("adapt","0");
    //
    addTorus.addArg("z","0");
    addTorus.addArg("r_major","0.7");
    addTorus.addArg("r_minor","0.2");
    addTorus.addArg("mu","0");
    addTorus.addArg("nZ","11");
    addTorus.addArg("adapt","0");
    //
    addWiggleSlab.addArg("z","0");
    addWiggleSlab.addArg("r","0.7");
    addWiggleSlab.addArg("end radius","0.05");
    addWiggleSlab.addArg("n wiggles","1");
    addWiggleSlab.addArg("amplitude","0.1");
    addWiggleSlab.addArg("mu","0");
    addWiggleSlab.addArg("nZ","11");
    addWiggleSlab.addArg("adapt","0");
    //
    OMsurfaces.addChoice(&setPhi,"n");
    OMsurfaces.addChoice(&addSlab,"s");
    OMsurfaces.addChoice(&addTube,"t");
    OMsurfaces.addChoice(&addBall,"b");
    OMsurfaces.addChoice(&addTorus,"l");
    OMsurfaces.addChoice(&addWiggleSlab,"w");
    OMsurfaces.addChoice(&punchHole,"hole",SELECTOR_HIDDEN);
    OMsurfaces.addChoice(&InputRequester::exitMenu,"x");
    
    doSolve.addArg("Saved solution file","none");
    
    qSurvey.addArg("N points","6");
    
    addSingular.addArg("State number","0");
    addSingular.addArg("Multiplier","50");
    
    setSingularEpsilon.addArg("epsilon","1e-5");
    
    punchHole.addArg("near x");
    punchHole.addArg("y");
    punchHole.addArg("z","0");
    punchHole.addArg("radius");
}

SystemConfiguration::~SystemConfiguration() {
    if(RSC) RSC->release();
    IncidentSource->release();
    TotalField->release();
    if(SS) delete SS;
}

void SystemConfiguration::initReactiveSet(unsigned int nPhi) {
    if(RSC) RSC->release();
    RSC = new MagRSCombiner(nPhi);
    RSC->retain();
    TotalField->addsource(RSC);
}


void SystemConfiguration::solve(const string& cfile) {
    assert(RSC); if(!RSC) return;
    if(SS) delete SS;
    SS = SymmetricSolver::cachedSolve(*RSC, cfile.size()?basedir+"/"+cfile:"");
}

void SystemConfiguration::calculate_result() {
    assert(SS); if(!SS) return;
    assert(RSC); if(!RSC) return;
    RSC->incidentState = RSC->getFullReactionTo(IncidentSource);
    SS->calculateResult(*RSC);
}

DVFunc1<2,double>* SystemConfiguration::adaptSurface(DVFunc1<2,double>* f, double pfixed, bool useTotal) const {
    assert(f);
    if(!pfixed || !IncidentSource->nSources()) return f;
    
    FieldEstimator2Dfrom3D fe(useTotal ? TotalField : IncidentSource);
    
    FieldAdaptiveSurface* FAS = new FieldAdaptiveSurface(*f);
    FAS->optimizeSpacing(fe,pfixed);
    return FAS;
}

void SystemConfiguration::measureFields(const string& xpath) const {
    // set up output files
    string fieldspath = basedir+"/"+xpath;
    makePath(fieldspath);
    std::ofstream fieldsout;
    std::ofstream statsout;
    std::ostream nullout(NULL);
    if(cell.saveGrid) fieldsout.open((fieldspath+"/Fieldmap.txt").c_str());
    statsout.open((fieldspath+"/Fieldstats.txt").c_str());
    
    // run analyzer
    if(Visualizable::vis_on) {
        TotalField->visualize();
        FA.visualizeSurvey(cell.ll,cell.ur,cell.vx,cell.vy,cell.vz);
    }
    FA.survey(cell.ll,cell.ur,cell.nx,cell.ny,cell.nz,statsout,cell.saveGrid?fieldsout:nullout);
    
    // cleanup
    if(cell.saveGrid) fieldsout.close();
    statsout.close();
}

void SystemConfiguration::writeInfo(const string& xpath) const {
    QFile qout;
    qout.insert("cell",cell.getInfo());
    qout.commit(basedir+"/"+xpath+"/GeomInfo.txt");
}

void SystemConfiguration::add_singular_state(unsigned int i, double c) {
    if(!SS || !RSC) return;
    SS->print_singular_values();
    VarVec<double> vi = SS->get_singular_values();
    assert(i<vi.size());
    VarVec<double> vs = SS->get_singular_vector(i);
    std::cout << "Adding state " << i << " with singular value " << vi[i]/vi[0] << " and magnitude " << vs.mag() << "\n";
    RSC->setFinalState(RSC->finalState + vs*c);
}


