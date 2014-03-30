#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *

sr = 0.05		# shield end radius
endgap = 0		# coil-to-main-shield end offset
sdr = 0.1		# shield-to-coil radial distance
plategap = 0.01	# main shield to endcaps offset

def make_asym_shortcoil_base(nm, r, dz=0, r0=0.70, gridt = 17, gride = 11):

	S = StudySetup(nm,r)
	S.nPhi = 32
	
	clen = 2.5
	crad = 1.0
	S.fields.append(CoilSpec(crad, clen, 15))	# r=1, l=2.5, N=15 cos theta coil
	S.fields[-1].dist = [-0.0025]
	S.fields[-1].ends = ["none","none"]			# "endless" wires
	
	srad = crad + sr + sdr
	shz = 0.5*clen + endgap
	S.shields.append(ShieldSpec(10000, gridt, [(-shz, srad),(shz, srad)], sr))
	
	S.shields.append(ShieldSpec(0, gride, [(-shz-sr-plategap, srad+sr)], sr))
	S.shields.append(ShieldSpec(0, gride, [(shz+sr+plategap,r0), (shz+sr+plategap, srad+sr)], sr))
	
	S.measGrid = (7,11,7)
	S.measCell = [(0.05,-0.20,-0.05+dz), (0.125,0.20,0.05+dz)]
	
	return S
	
def make_bothends_closed(nm, r, dz=0, gridt = 25, gride = 17):
	"""Both ends closed superconductor"""
	S = make_asym_shortcoil_base(nm,r,dz, gridt=gridt, gride=gride)
	S.shields[-1] = ShieldSpec(0, gride, [(-S.shields[-2].ends[0][0],S.shields[-2].ends[0][1])], S.shields[-2].rthick)
	return S

def make_smallholes_annular(nm,r,dz=0):
	"""More 'realistic' endcap with small central hole and annular ring"""
	S = make_asym_shortcoil_base(nm,r,dz,r0=0.60)
	
	shz = S.shields[-1].ends[0][0]
	S.shields.append(ShieldSpec(0, 17, [(shz,0.0635),(shz,0.52)], sr))
	
	return S

#make_setup = make_asym_shortcoil_base
#stname = "ShortCoil"

make_setup = make_bothends_closed
stname = "ClosedCoil"

#make_setup = make_smallholes_annular
#stname = "SC_SmallAnnular"


def DistortionScan():
	SS = StudyScan()
	for (n,r) in enumerate(unifrange(-0.01, 0.01, 8, True)):
		S = make_setup(stname+"_Distortion",r)
		S.solfl = "../MC"
		S.fields[-1].dist = [r]
		if False: #n==0:
			os.system(S.make_cmd("RotationShield_Vis"))
		else:
			SS.fsimlist.write(S.make_cmd())
		print S.make_cmd();
	SS.run()

def MovingCellScan():
	SS = StudyScan()
	for (n,r) in enumerate(unifrange(-1, 1, 8, True)):
		S = make_setup(stname+"_MovingCell",r,dz=r)
		S.solfl = "../MC"
		if n == 0:
			os.system(S.make_cmd("RotationShield_Vis"))
		else:
			SS.fsimlist.write(S.make_cmd())
	SS.run()

def GriddingScan():
	SS = StudyScan()
	for r in range(0,8):
		S = make_setup(stname+"_Gridding", r, gridt=11+2*r, gride=7+r)
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def PhiGriddingScan():
	SS = StudyScan()
	for r in range(0,17)[::-1]:
		nPhi = int(2**(3+r/4.))
		S = make_setup(stname+"_PhiGridding", nPhi, gridt=29, gride=16)
		S.nPhi = nPhi
		SS.fsimlist.write(S.make_cmd())
	SS.run()


#
# want average gradients < 0.1 uG/cm
#

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
		PhiGriddingScan()

	outdir = os.environ["ROTSHIELD_OUT"]

	VPP = None

	if options.plot:
	
		#VPP = VarParamPlotter(outdir+"/"+stname+"_MovingCell")
		#VPP.keypos = "tr"
		#VPP.setupGraph("Cell center z [m]")
		
		#VPP = VarParamPlotter(outdir+"/"+stname+"_Distortion")
		#VPP.keypos = "tl"
		#VPP.setupGraph("Distortion parameter $a$")
	
	
		#VPP = VarParamPlotter(outdir+"/"+stname+"_ECdz")
		#VPP.keypos = "tl"
		#VPP.setupGraph("Wires to endcap distance [m]")

		#VPP = VarParamPlotter(outdir+"/"+stname+"_Gridding")
		#VPP.keypos = "tr"
		#VPP.setupGraph("longitudinal grid divisions $N_Z$ ($N_\\phi = 32$)")
		
		VPP = VarParamPlotter(outdir+"/"+stname+"_PhiGridding")
		VPP.keypos = "tr"
		VPP.setupGraph("radial grid divisions $N_\\phi$ ($N_Z = 61$)", logx=True, xrange=(7,150))
		

		if VPP is not None:
			VPP.makePlot(PGlist=[PG_n(),PG_3He()])
			VPP.outputPlot()
			

	if options.fplt:
		if VPP is not None:
			for r,FI in VPP.datlist:
				FI.BC.plotFields(2,0,1,2)
		else:
			
			FI = FieldInfo(outdir+"/"+stname+"_ECdz/X_0.050000/")
			#FI = FieldInfo(outdir+"/foo/")
			
			if 1:
				FI.BC.plotFields(2,0,1,2) # Bz along z
				FI.BC.plotFields(0,0,2,1) # Bx along y
				FI.BC.plotFields(0,0,1,2) # Bx along z
			if 1:
				FI.BC.plotCellProjection(1,0)
				FI.BC.plotCellProjection(1,2)
				FI.BC.plotCellProjection(0,2)


