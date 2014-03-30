/* 
 * CosThetaBuilder.hh, part of the RotationShield program
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

#ifndef COSTHETABUILDER_HH
#define COSTHETABUILDER_HH 1

#include "MixedSource.hh"
#include "Typedefs.hh"
#include "QFile.hh"
#include "ControlMenu.hh"
#include <vector>

/// base class for positioning wires
class AnglePositioner {
public:
	/// constructor
	AnglePositioner() {}
	/// destructor
	virtual ~AnglePositioner() {}
	/// get wire angle
	virtual double angle(unsigned int i, unsigned int ncoils) const = 0;
	/// get information in Stringmap form
	virtual Stringmap getInfo() const = 0;
};

/// position cos theta using fourier-series offsets
class ShiftPositioner: public AnglePositioner {
public:
	/// constructor
	ShiftPositioner(const VarVec<double>& v = VarVec<double>(0)): shift(v) {}
	/// get wire angle
	virtual double angle(unsigned int i, unsigned int ncoils) const;
	/// get information in Stringmap form
	virtual Stringmap getInfo() const;
	
	VarVec<double> shift; //< shifting parameters
};

//--------------------------------------------------------------------

/// base class for translating wires
class EndTranslator {
public:
	/// constructor
	EndTranslator() {}
	/// destructor
	virtual ~EndTranslator() {}
	/// get translation
	virtual vec3 trans(int N, int i, int a1, int a2, int a3) = 0;
	/// get information in Stringmap form
	virtual Stringmap getInfo() const { return Stringmap(); }
};

/// simple vector translation
class VecTrans: public EndTranslator {
public:
	/// constructor
	VecTrans(vec3 t = vec3()): trans1(t), trans2(t) {}
	/// get translation
	virtual vec3 trans(int, int, int a1, int, int) { return a1>0?trans1:trans2; }
	/// get information in Stringmap form
	virtual Stringmap getInfo() const;
	
	vec3 trans1;	//< translation for one end
	vec3 trans2;	//< translation for other end
};

//--------------------------------------------------------------------

class CosThetaBuilder: public StreamInteractor {
public:
	/// constructor
	CosThetaBuilder(unsigned int n = 15, double r = 0.5, double l = 1.0,
					AnglePositioner* ap = new ShiftPositioner(), EndTranslator* et = NULL);
	
	/// build coil into MixedSource
	void buildCoil(MixedSource& M);
	/// write info to QFile
	void writeInfo(QFile& qOut) const;
	
	unsigned int ncoils;		//< number of saddle coils per half
	float radius;				//< coil radius
	float length;				//< coil length
	
	AnglePositioner* AP;
	EndTranslator* ET;
	
	/// get specified wire end position
	vec3 getEndp(unsigned int n, bool xside, bool yside, bool zside);
	
	enum capType {
		CAP_LINE,		//< line segments between endpoints
		CAP_STRAIGHT,	//< rectangular coils straight across ends
		CAP_ARC,		//< curved arcs between endpoints
		CAP_NONE		//< no endcap wires constructed
	};
	capType myCap[2];	//< how to construct endcaps on each side
	unsigned int nArc;
	
	// menu user interface
	InputRequester setGeometry;
	InputRequester setDistortion;
	NameSelector selectEndcap;
	InputRequester setEndcaps;
	OptionsMenu OMcoil;
	
	
protected:
	
	std::vector<vec3> endp;
	
	void buildEndpoints();
	void buildSides(MixedSource& M);
		
	void buildArcCap(MixedSource& M, unsigned int zside, unsigned int nseg=1);
	void buildLineCap(MixedSource& M, unsigned int zside);
	void buildStraightCap(MixedSource& M, unsigned int zside);
	
	void buildMixedCaps(MixedSource& M, float rinner);
};

#endif
