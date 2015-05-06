#!/usr/bin/python

from ShieldStudyLauncher import *

def make_aCORN(nm,r):
    S = StudySetup(nm,r)
    S.nPhi = 32
    
    ncoils = 25
    crad = 0.3
    dcoils = 0.1
    cj = 1.0
    coilj = (lambda i: (i < ncoils-3)*1. + (i == ncoils-3)*1.7)
    coilpos = [ (coilj(i), (i-0.5*ncoils+0.5)*dcoils) for i in range(ncoils)]
    #if ncoils%2:
    #    coilpos.append((0.5*cj, 0));
    for z in coilpos:
        print "Coil at", z[1], "I =", z[0]
        S.fields.append(RingSpec(crad, z[1], z[0]))
        
    #S.fields.append(commandSpec("s 1"));

    gride = 12
    sdz = 0.07
    shz = 0.5*ncoils*dcoils + sdz
    s_ir = 0.2
    s_or = 0.75
    mu = 1000
    S.shields.append(ShieldSpec(mu, gride, [(shz, s_or)], sdz, 0.7))
    S.shields.append(ShieldSpec(mu, gride, [(-shz, s_or)], sdz, 0.7))
    #S.shields.append(ShieldSpec(mu, 1.5*gride, [(-shz+2*sdz, s_or), (shz-2*sdz, s_or)], sdz, 0.7))
    mid_z = (0.5*ncoils-2)*dcoils
    S.shields.append(ShieldSpec(mu, 1.5*gride, [(mid_z, crad), (mid_z, s_or)], 0.04, 0.7))
    
    
    S.measCell = [(0,0,0),(0,0,1.5)]
    S.measGrid = (1,1,151)
    return S

if __name__ == "__main__":
    S = make_aCORN("aCORN",0)
    print S.make_cmd()
    