/* 
 * Studies.hh, part of the RotationShield program
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

/// \file Studies.hh \brief Utility classes for setting up Cos Theta Coil studies

#ifndef STUDIES_HH
/// Make sure this header is only loaded once
#define STUDIES_HH

#include "MixedSource.hh"
#include "SymmetricSolver.hh"
#include "CosThetaBuilder.hh"
#include "FieldAnalyzer.hh"
#include "MagRS.hh"
#include "ControlMenu.hh"
#include "QFile.hh"
#include "DVFunc.hh"
#include <string>

/// Data-sampling cell specification
class fieldCell {
public:
    /// default constructor
    fieldCell(): ll(-.2,-.2,-.2), ur(.2,.2,.2), nx(5), ny(5), nz(5), vx(2), vy(3), vz(2), saveGrid(true) {}
    /// cell info as Stringmap
    Stringmap getInfo() const;
    
    vec3 ll;                    ///< lower corner
    vec3 ur;                    ///< upper corner
    unsigned int nx,ny,nz;      ///< number of points to sample along each axis
    unsigned int vx,vy,vz;      ///< number of points to visualize along each axis
    bool saveGrid;              ///< whether to save gridded output on survey
};

/// Utility class for configuring field measurement
class SystemConfiguration: public StreamInteractor {
public:
    /// constructor
    SystemConfiguration();
    /// destructor
    ~SystemConfiguration();
    
    /// initialize symmetric boundary system
    void initReactiveSet(unsigned int nPhi);
    
    /// adapt a geometry to current field source
    DVFunc1<2,double>* adaptSurface(DVFunc1<2,double>* f, double pfixed, bool useTotal = false) const;
    /// solve surface interaction system
    void solve(const string& cfile = "");
    /// calculate applied response
    void calculate_result();
    
    /// measure fields, writing to given output file
    void measureFields(const string& xpath="") const;
    /// write measurement info to file at basedir/xpath/GeomInfo.txt
    void writeInfo(const string& xpath="") const;
    
    /// display singular solution
    void add_singular_state(unsigned int i, double c);
    
    
    string basedir;                ///< base directory for IO operations
    MagRSCombiner* RSC;                 ///< reacting boundary condition surfaces
    MagRSCombiner* PTB;                 ///< perturbation reaction
    MagExtField* IncidentSource;        ///< incident field source
    MixedSource* TotalField;            ///< incident + reacting field
    FieldAnalyzer FA;
    fieldCell cell;                     ///< field measurement cell
    SymmetricSolver* SS;                ///< system solver
    CosThetaBuilder CTB;                ///< cos theta field coil builder
    bool updateVis = true;              ///< whether to update visualization
    
    
    //------------------------------------
    // menu-driven user interface elements
    
    InputRequester outDir;
    
    InputRequester setFCrange;
    InputRequester setFCgrid;
    InputRequester setSaveGrid;
    OptionsMenu OMcell;
    
    InputRequester addLineCurrent;
    InputRequester addUnifB;
    InputRequester clearIncident;
    InputRequester buildCosThetaExit;
    InputRequester addRingCurrent;
    InputRequester symmetrizeField;
    InputRequester toggleUpdateVis;
    OptionsMenu OMfieldsrc;
    
    InputRequester setPhi;
    InputRequester addSlab;
    InputRequester addTube;
    InputRequester addBall;
    InputRequester addTorus;
    InputRequester addWiggleSlab;
    OptionsMenu OMsurfaces;
    
    InputRequester doSolve;
    InputRequester doApply;
    InputRequester zeroResponse;
    InputRequester doMeas;
    InputRequester qSurvey;
    InputRequester addSingular;
    InputRequester setSingularEpsilon;

    // hidden advanced functions
    InputRequester punchHole;
    InputRequester equilibratePtb;
};

#endif
