#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *

endgap = 0        # coil-to-main-shield end offset
sdr = 0.1        # shield-to-coil radial distance

def make_asym_shortcoil_base(nm, r,
                             dz=-0.75,          # cell z offset
                             r0=0.50,           # top hole radius
                             gridt = 58,        # tube Nz grid
                             gride = 32,        # endcap Nz grid
                             gridhx = 0,        # extra gridding on hole side
                             plategap = 0.01,   # shield-to-cap gap
                             coilgap = 0.01,    # cap-to-coil-end gap
                             sr = 0.02          # volume end radius
                             ):

    S = StudySetup(nm,r)
    S.nPhi = 64
    
    clen = 2.5
    crad = 1.0
    S.fields.append(CoilSpec(crad, clen + 2*(plategap-coilgap), 15))    # r=1, l=2.5, N=15 cos theta coil
    S.fields[-1].dist = [-0.0032]                # Closed: -.0032; smallhole -0.0033; r=70cm -.0005; flattener -0.0036
    S.fields[-1].ends = ["none","none"]            # "endless" wires
    
    srad = crad + sr + sdr
    shz = 0.5*clen + endgap
    S.shields.append(ShieldSpec(10000, gridt, [(-shz, srad),(shz, srad)], sr))
    
    S.shields.append(ShieldSpec(0, gride, [[-shz-sr-plategap, srad+sr]], sr))
    S.shields.append(ShieldSpec(0, gride+gridhx, [(shz+sr+plategap,r0), (shz+sr+plategap, srad+sr)], sr))
    
    S.measGrid = (7,11,7)
    S.measCell = [(0.05,-0.20,-0.05+dz), (0.125,0.20,0.05+dz)]
    
    return S
    
def make_bothends_closed(nm, r, **kwargs):
    """Both ends closed superconductor"""
    S = make_asym_shortcoil_base(nm, r, **kwargs)
    S.shields[-1] = ShieldSpec(0, S.shields[-2].nseg, [(-S.shields[-2].ends[0][0],S.shields[-2].ends[0][1])], S.shields[-2].rthick)
    return S

def make_smallholes_annular(nm, r, **kwargs):
    """More 'realistic' endcap with small central hole and annular ring"""
    S = make_asym_shortcoil_base(nm, r, r0=0.60, **kwargs)
    
    shz = S.shields[-1].ends[0][0]
    S.shields.append(ShieldSpec(0, S.shields[-1].nseg, [(shz,0.0635),(shz,0.52)], S.shields[-1].rthick))
    
    return S

def make_holey(nm,r,**kwargs):
    """'Realistic' endcap with many hole perturbations"""
    S = make_bothends_closed(nm, r, **kwargs)
    
    S.shields[-1].holes.append(ShieldHoleSpec((0,0,0), .127*0.5))
    S.shields[-1].holes.append(ShieldHoleSpec((.2963, -.2963, 0), .127*0.5))
    
    S.shields[-1].holes.append(ShieldHoleSpec((-.4572, 0, 0), .1778*0.5))
    S.shields[-1].holes.append(ShieldHoleSpec((.2540, 0, 0), .1778*0.5))

    r0 = 1.1176*0.5
    a = .0762*0.5
    for i in range(5):
        th = 2*pi*i/36.
        for sgx in (-1,1):
            for sgy in (-1,1):
                S.shields[-1].holes.append(ShieldHoleSpec((r0*cos(0.61+th)*sgx, r0*sin(0.61+th)*sgy, 0), a))
    S.shields[-1].holes.append(ShieldHoleSpec((r0, 0, 0), a))
    S.shields[-1].holes.append(ShieldHoleSpec((-0.0191, .5482, 0), a))
    S.shields[-1].holes.append(ShieldHoleSpec((-0.0191, -.5482, 0), a))
        
    r0 = 1.32*0.5
    for i in range(6):
        th = 2*pi*i/6.
        S.shields[-1].holes.append(ShieldHoleSpec((r0*cos(th), r0*sin(th), 0), a))
    
    return S

def make_wiggley(nm, r, nwiggles = 3, wsize = 0.02, **kwargs):
    S = make_bothends_closed(nm, r, **kwargs)
    S.shields[-2].nwiggles = nwiggles
    S.shields[-2].wigglesize = wsize
    if nwiggles%2:
        S.shields[-2].ends[0][0] += 2*wsize
    return S

def make_outercoil(nm,r,**kwargs):
    """Big-hole shield with outside cos theta flattener"""
    S = make_asym_shortcoil_base(nm, r, **kwargs)
    S.fields.append(CoilSpec(1.0, 0.5, 15))
    S.fields[-1].offset = (0,0, 0.25 + 1.25 + 2*0.01 + 0.04)
    S.fields[-1].ends = ["line","none"]
    S.fields[-1].j = 1.29

    return S;

#make_setup = make_outercoil
#stname = "ShortCoil_50cm"
#stname = "SC_SmallHole"

make_setup = make_bothends_closed
stname = "ClosedCoil"

#make_setup = make_smallholes_annular
#stname = "SC_SmallAnnular"

#make_setup = make_holey
#stname = "SC_manyholes"

#make_setup = make_wiggley
#stname = "SC_Wiggles"

def DistortionScan():
    SS = StudyScan()
    for (n,r) in enumerate(unifrange(-0.02, 0, 8, True)):
        S = make_setup(stname+"/Distortion",r)
        S.solfl = "../../SC"
        for f in S.fields:
            f.dist = [r, -0.020, -0.030]
        if n==-1:
            os.system(S.make_cmd("RotationShield_Vis"))
        else:
            SS.fsimlist.write(S.make_cmd())
    SS.run()

def MovingCellScan():
    SS = StudyScan()
    for (n,r) in enumerate(unifrange(-1, 1, 32, True)):
        S = make_setup(stname+"/MovingCell",r,dz=r)
        S.solfl = "../../SC"
        if n == 0:
            os.system(S.make_cmd("RotationShield_Vis"))
        else:
            SS.fsimlist.write(S.make_cmd())
    SS.run()

def GriddingScan():
    SS = StudyScan()
    for r in range(0,31)[::-1]:
        S = make_setup(stname+"/Gridding", r, gridt=11+2*r, gride=7+r)
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def HoleGridScan():
    SS = StudyScan()
    for r in range(0,24)[::3][::-1]:
        S = make_setup(stname+"/HoleGridding", r, gridhx = r)
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def PhiGriddingScan():
    SS = StudyScan()
    for r in range(0,17)[::-1]:
        nPhi = int(2**(3+r/4.))
        S = make_setup(stname+"/PhiGridding", nPhi, gridt=29, gride=16)
        S.nPhi = nPhi
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def PlateGapScan():
    SS = StudyScan()
    for r in unifrange(.01, .1, 8, True) + unifrange(.001, .015, 8, True):
        S = make_setup(stname+"/PlateGap", r, dz=-1, plategap = r)
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def EdgeRadiusScan():
    SS = StudyScan()
    for r in  unifrange(.01, .1, 8) + unifrange(.001, .015, 8):
        S = make_setup(stname+"/EdgeRadius", r, dz=0, sr = r)
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def HoleSizeScan():
    SS = StudyScan()
    for r in  unifrange(.02, .5, 8):
        S = make_setup(stname+"/HoleSize", r, r0=r)
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def HolePertScan():
    SS = StudyScan()
    for (n,r) in  enumerate(unifrange(.11, .5, 8)): #unifrange(.02, .5, 8):
        S = make_setup(stname+"/HolePert", r)
        S.shields[-1].holes.append(ShieldHoleSpec((0,0,0),r))
        S.solfl = "../SC"
        if n == -1:
            os.system(S.make_cmd("RotationShield_Vis"))
        else:
            SS.fsimlist.write(S.make_cmd())
    SS.run()

def RotateEnd():
    SS = StudyScan()
    for r in  unifrange(0, 2*pi, 32):
        S = make_setup(stname+"/EndRot", r)
        S.shields[-1].rot_pert(r)
        S.solfl = "../../SC"
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def Outer_Jscan():
    SS = StudyScan()
    for r in  unifrange(0, 2, 32):
        S = make_setup(stname+"/OuterJ", r)
        S.fields[-1].j = r
        S.solfl = "../../SC"
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def WiggleScan():
    SS = StudyScan()
    nw = 3
    for r in  unifrange(-0.02, 0.02, 7):
        S = make_setup(stname+"/WiggleSize", r, wsize = r, nwiggles=nw, dz = -0.5)
        #S.solfl = "SC_%i"%nw
        SS.fsimlist.write(S.make_cmd())
    SS.run()

def CoilEndScan():
    SS = StudyScan()
    for (n,r) in  enumerate(unifrange(0.001, 0.02, 8)):
        S = make_setup(stname+"/WireDistance", r, dz=-0.5, gride = 64, coilgap = r)
        S.fields[-1].ends = None
        S.shields[-1].p = S.shields[-2].p = 0.5
        S.solfl = "SC"
        #if n == 7:
        #    os.system(S.make_cmd("RotationShield_Vis"))
        #else:
        SS.fsimlist.write(S.make_cmd())
    SS.run()



if __name__=="__main__":

    parser = OptionParser()
    parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
    parser.add_option("--scan", dest="scan", action="store_true", default=False, help="run simulation jobs")
    parser.add_option("--plot", dest="plot", action="store_true", default=False, help="plot simulation results")
    parser.add_option("--fplt", dest="fplt", action="store_true", default=False, help="plot field shape")
    
    options, args = parser.parse_args()
    if options.kill:
        os.system("killall -9 parallel")
        os.system("killall -9 RotationShield")
        exit(0)
    
    if options.scan:
        
        #MovingCellScan()
        #DistortionScan()
        
        #GriddingScan()
        #PhiGriddingScan()
        #PlateGapScan()
        #EdgeRadiusScan()
        #HoleSizeScan()
        #HolePertScan()
        #HoleGridScan()
        #RotateEnd()
        #Outer_Jscan()
        #WiggleScan()
        CoilEndScan()

    outdir = os.environ["ROTSHIELD_OUT"]

    VPP = None

    if options.plot:
    
        #VPP = VarParamPlotter(outdir+"/"+stname+"/MovingCell")
        #VPP.keypos = "tc"
        #VPP.setupGraph("Cell center z [m]")
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/Distortion")
        #VPP.keypos = "tl"
        #VPP.setupGraph("Distortion parameter $a_1, -0.02, -0.03$")
    
    
        #VPP = VarParamPlotter(outdir+"/"+stname+"/ECdz")
        #VPP.keypos = "tl"
        #VPP.setupGraph("Wires to endcap distance [m]")

        #VPP = VarParamPlotter(outdir+"/"+stname+"/Gridding")
        #VPP.keypos = "tr"
        #VPP.setupGraph("longitudinal grid divisions $N_Z$ ($N_\\phi = 32$)", logy2=True, xtrans=(lambda x: 25+4*x))
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/PhiGridding")
        #VPP.keypos = "tr"
        #VPP.setupGraph("radial grid divisions $N_\\phi$ ($N_Z = 61$)", logx=True, xrange=(7,150))
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/PlateGap")
        #VPP.keypos = "tl"
        #VPP.setupGraph("shield-to-endplate gap [cm]", xtrans=(lambda x: 100*x))
    
        #VPP = VarParamPlotter(outdir+"/"+stname+"/EdgeRadius")
        #VPP.keypos = "tc"
        #VPP.setupGraph("Solid edge radius [cm]", xtrans=(lambda x: 100*x))

        #VPP = VarParamPlotter(outdir+"/"+stname+"/HoleSize")
        #VPP.keypos = "tl"
        #VPP.setupGraph("End hole radius [cm]", xtrans=(lambda x: 100*x), yrange=(-8,6), y2range=(0,15))
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/HolePert")
        #VPP.keypos = "tl"
        #VPP.setupGraph("perturbation hole radius [cm]", xtrans=(lambda x: 100*x), yrange=(-8,6), y2range=(0,15))
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/HoleGridding")
        #VPP.keypos = "tr"
        #VPP.setupGraph("additional grid divisions")
        
        #VPP = VarParamPlotter(outdir+"/"+stname+"/EndRot")
        #VPP.keypos = "tc"
        #VPP.setupGraph("Endcap rotation angle [degrees]", xrange=(0,360), xtrans=(lambda x: 180*x/pi))
    
        #VPP = VarParamPlotter(outdir+"/"+stname+"/OuterJ")
        #VPP.keypos = "tc"
        #VPP.setupGraph("External coil relative current")

        #VPP = VarParamPlotter(outdir+"/"+stname+"/WiggleSize")
        #VPP.keypos = "tc"
        #VPP.setupGraph("Bottom plate wiggle amplitude [cm]", xtrans=(lambda x: 100*x))

        VPP = VarParamPlotter(outdir+"/"+stname+"/WireDistance")
        VPP.keypos = "tc"
        VPP.setupGraph("Endcap-to-coil gap [cm]", xtrans=(lambda x: 100*x))

        if VPP is not None:
            VPP.makePlot(PGlist=[PG_n(),PG_3He()])
            VPP.outputPlot()
            

    if options.fplt:
        if VPP is not None:
            for r,FI in VPP.datlist:
                FI.BC.plotFields(2,0,1,2)
        else:
            
            FI = FieldInfo(outdir+"/"+stname+"/ECdz/X_0.050000/")
            #FI = FieldInfo(outdir+"/foo/")
            
            if 1:
                FI.BC.plotFields(2,0,1,2) # Bz along z
                FI.BC.plotFields(0,0,2,1) # Bx along y
                FI.BC.plotFields(0,0,1,2) # Bx along z
            if 1:
                FI.BC.plotCellProjection(1,0)
                FI.BC.plotCellProjection(1,2)
                FI.BC.plotCellProjection(0,2)


