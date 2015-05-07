#!/usr/bin/python

from ShieldStudyLauncher import *

def make_aCORN(nm,r):
    S = StudySetup(nm,r)
    S.nPhi = 32
    
    ncoils = 25                 # number of coils
    nsimwind = 3                # number of simulated winding radii
    nsimz = 1                   # number of simulated winding heights
    crad0 = 0.22225             # main coil inner radius
    crad1 = 0.385               # main coil outer radius
    dcoils = 0.12               # spacing between coils
    cdz = 0.02                  # height of each coil
    cj = 9.5*121/1000;          # 9.5A in 121 turns current; scaled down by 1000
    #coilj = (lambda i: (i < ncoils-3)*1. + (i == ncoils-3)*1.7)
    coilpos = [ (cj, (i-0.5*ncoils+0.5)*dcoils) for i in range(ncoils)]
    for z in coilpos:
        print "Coil at", z[1], "I =", z[0]
        for i in range(nsimwind):
            crad = crad0 + (i+0.5)/nsimwind*(crad1 - crad0)
            for iz in range(nsimz):
                cz = z[1] + (iz + 0.5)/nsimz*cdz
                S.fields.append(RingSpec(crad, cz, z[0]/(nsimz*nsimwind)))

    gride = 17          # flux return gridding
    sdz = 0.0254/2      # flux return thickness (really 0.0254/2)
    shz = 0.5*ncoils*dcoils + sdz
    s_ir = 0.088        # flux return inner radius
    s_or = 0.80         # flux return outer radius
    mu = 1000           # flux return mu
    frb_z = coilpos[0][1]-0.0581-sdz    # bottom flux return z
    frt_z = coilpos[-1][1]+0.0639+sdz   # top flux return z
    
    ########### Bottom
    #S.shields.append(ShieldSpec(mu, gride, [(frb_z, s_or)], sdz, 0.7))
    S.shields.append(ShieldSpec(mu, 1.5*gride, [(frb_z, s_ir), (frb_z, s_or)], sdz, 0.7))
    
    ########### Top
    S.shields.append(ShieldSpec(mu, gride, [(frt_z, s_or)], sdz, 0.7))
    
    ########### Flux return side can
    #S.shields.append(ShieldSpec(mu, 1.5*gride, [(frb_z+1.5*sdz, s_or), (frt_z-1.5*sdz, s_or)], sdz, 0.7))
    
    #mid_z = (0.5*ncoils-2)*dcoils
    #S.shields.append(ShieldSpec(mu, 1.5*gride, [(mid_z, crad), (mid_z, s_or)], 0.04, 0.7))
    
    
    S.measCell = [(0,0,0),(0,0,2.0)]
    S.measGrid = (1,1,101)
    return S

if __name__ == "__main__":
    S = make_aCORN("aCORN",0)
    print S.make_cmd()
    