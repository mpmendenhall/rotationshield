#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *


def make_asym_shortcoil_base(nm,r,dz=-0.783,r0=0.70):

	S = StudySetup(nm,r)
	
	clen = 2.5
	crad = 1.0
	S.fields.append(CoilSpec(crad, clen, 15))	# r=1, l=2.5, N=15 cos theta coil
	S.fields[-1].dist = [-0.0096]

	sr = 0.05		# shield end radius
	endgap = 0.02	# coil-to-main-shield end offset
	sdr = 0.1		# shield-to-coil distance
	srad = crad + sr + sdr
	shz = 0.5*clen + endgap
	S.shields.append(ShieldSpec(10000, 17, [(-shz, srad),(shz, srad)], sr))
			
	S.measGrid = (7,11,7)
	S.measCell = [(0.05,-0.20,-0.05+dz), (0.125,0.20,0.05+dz)]
	
	return S
	

def make_smallholes_annular(nm,r,dz=0):
	"""More 'realistic' endcap with small central hole and annular ring"""
	S = make_asym_shortcoil_base(nm,r,dz,r0=0.60)
	S.set_csgeom(clen = 2.5, crad = 1.0, dlen = 0.01, drad = 0.1)
	
	S.shsects.append(ShieldSection(0,6,12,"pax","pax"))
	S.shsects[-1].off = ([0,0.52],[0,0.0635])
	
	S.shsects[0].off[0][1] = 0.50	# bottom end hole
	#S.shsects[0].cseg = 24
	#S.shsects[0].vseg = 12
		
	S.cend = ["none","none"]
	S.dist = [-0.0184]
	
	return S

make_setup = make_asym_shortcoil_base
stname = "ShortCoil"





def ECtoWiresScan():
	SS = StudyScan()
	for r in unifrange(0.001, .05, 7, True):
		S = make_setup(stname+"_ECdz",r)
		S.set_csgeom(clen = 2.5, crad = 1.0, dlen = r, drad = 0.1)
		#print S.make_cmd();
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def DistortionScan():
	SS = StudyScan()
	for (n,r) in enumerate(unifrange(-0.04, 0, 8, True)):
		S = make_setup(stname+"_Distortion",r)
		S.solfl = "../SCD"
		S.dist = [r]
		if not n:
			os.system(S.make_cmd("RotationShield_Vis"))
		else:
			SS.fsimlist.write(S.make_cmd())
		print S.make_cmd();
	SS.run()

def MovingCellScan():
	SS = StudyScan()
	for (n,r) in enumerate(unifrange(-.5, 0.5, 8, True)):
		S = make_setup(stname+"_MovingCell",r,dz=r)
		S.solfl = "../MC"
		if not n:
			os.system(S.make_cmd("RotationShield_Vis"))
		else:
			SS.fsimlist.write(S.make_cmd())
	SS.run()


def PosBlockScan():
	SS = StudyScan()
	for r in unifrange(0, 0.7, 7, True):
		S = make_setup(stname+"_PosBlock",r)
		#S.shsects[1].off[1][0] = S.shsects[2].off[0][0] = S.shsects[2].off[1][0] = r
		S.shsects[3].off[0][1] = r
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def GriddingScan():
	SS = StudyScan()
	for r in unifrange(12, 27, 16, True):
		S = make_setup(stname+"stname_Gridding",r)
		#S.shsects[0].cseg = S.shsects[2].cseg = r
		S.shsects[1].cseg = r
		S.shsects[1].vseg = 0
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
		DistortionScan()
		
		#ECtoWiresScan()
				
		#NegApertureScan()
		
		#GriddingScan()

	outdir = os.environ["ROTSHIELD_OUT"]

	VPP = None

	if options.plot:
	
		#VPP = VarParamPlotter(outdir+"/"+stname+"_MovingCell")
		#VPP.keypos = "tc"
		#VPP.setupGraph("Cell center z [m]")
		
		VPP = VarParamPlotter(outdir+"/"+stname+"_Distortion")
		VPP.keypos = "tr"
		VPP.setupGraph("Distortion parameter $a$")
	
	
		#VPP = VarParamPlotter(outdir+"/"+stname+"_ECdz")
		#VPP.keypos = "tl"
		#VPP.setupGraph("Wires to endcap distance [m]")

		#VPP = VarParamPlotter(outdir+"/"+stname+"_NegAperture")
		#VPP.keypos = "tc"
		#VPP.setupGraph("Baseplate aperture [m]")

		#VPP = VarParamPlotter(outdir+"/"+stname+"_Gridding")
		#VPP.setupGraph("Endcap grid segments")

		#VPP = VarParamPlotter(outdir+"/"+stname+"_Cone")
		#VPP.setupGraph("Endcap cone offset [m]")

		#VPP = VarParamPlotter(outdir+"/"+stname+"_PosBlock")
		#VPP.setupGraph("Inner disc z offset [m]")
		
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


