#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *


def make_asym_shortcoil_base(nm,r,dz=-0.819):

	S = StudySetup(nm,r)
	S.set_csgeom(clen = 2.5, crad = 1.0, dlen = 0.05, drad = 0.1)
	
	S.shsects.append(ShieldSection(0,6,12,"nax","nrd"))
	S.shsects.append(ShieldSection(10000,10,20,"nrd","prd"))
	S.shsects.append(ShieldSection(0,6,12,"prd","pax"))
	S.shsects[2].off[1][1] = 0.70
		
	S.measGrid = (7,11,7)
	S.measCell = [(0.05,-0.20,-0.05+dz), (0.125,0.20,0.05+dz)]
	S.dist = [-0.0116]
		
	return S

def make_asym_shortcoil_throughend(nm,r,dz=-0.10): #2751):

	S = make_asym_shortcoil_base(nm,r,dz)
	S.set_csgeom(clen = 2.5, crad = 1.0, dlen = 0.01, drad = 0.1)
	
	S.shsects[2].off[1][1] = 0.95	# just cover endcap
	# move close to endcap
	S.shsects[1].off[1][0] = S.shsects[2].off[0][0] = S.shsects[2].off[1][0] = -0.005
	
	# inner plug
	S.shsects.append(ShieldSection(0,6,12,"pax","pax"))
	S.shsects[3].off[0][1] = 0.65
	
	#S.shsects[2].off[1][1] = 0
	#S.shsects.append(ShieldSection(0,6,12,"prd","pax"))
	#S.shsects[3].off[1][1] = 0.95
	
	S.cend = ["none",None]
	S.dist = [-0.0080]
	
	return S



make_setup = make_asym_shortcoil_throughend

def BothApertureScan():
	SS = StudyScan()
	for r in unifrange(0.05, 1, 7, True):
		S = make_setup("ShortCoil_Aperture",r)
		S.ecapR = [r,r]
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def ECtoWiresScan():
	SS = StudyScan()
	for r in unifrange(0.001, .05, 7, True):
		S = make_setup("ShortCoil_ECdz",r)
		S.set_csgeom(clen = 2.5, crad = 1.0, dlen = r, drad = 0.1)
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def DistortionScan():
	SS = StudyScan()
	for r in unifrange(-0.02, 0, 7, True):
		S = make_setup("ShortCoil_Distortion",r)
		S.dist = [r]
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def MovingCellScan():
	SS = StudyScan()
	for r in unifrange(-1, 0.5, 7, True):
		S = make_setup("ShortCoil_MovingCell",r,dz=r)
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def LengthScan():
	SS = StudyScan()
	for r in unifrange(2.5, 10, 7, True):
		S = make_setup("ShortCoil_Length",r)
		S.set_csgeom(clen = r, crad = 1.0, dlen = 0.01, drad = 0.1)
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def NegApertureScan():
	SS = StudyScan()
	for r in unifrange(0, 0.2, 7, True):
		S = make_setup("ShortCoil_NegAperture",r)
		S.shsects[0].off[0][1] = r
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def PosBlockScan():
	SS = StudyScan()
	for r in unifrange(0, 0.7, 7, True):
		S = make_setup("ShortCoil_PosBlock",r)
		#S.shsects[1].off[1][0] = S.shsects[2].off[0][0] = S.shsects[2].off[1][0] = r
		S.shsects[3].off[0][1] = r
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def GriddingScan():
	SS = StudyScan()
	for r in unifrange(12, 27, 16, True):
		S = make_setup("ShortCoil_Gridding",r)
		#S.shsects[0].cseg = S.shsects[2].cseg = r
		S.shsects[1].cseg = r
		S.shsects[1].vseg = 0
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def SuperAScan():
	SS = StudyScan()
	for r in unifrange(-0.02, 0, 7, True):
		S = make_setup("ShortCoil_a5",r)
		S.dist = [-0.028,0,-0.3,0,r]
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def ConeScan():
	SS = StudyScan()
	for r in unifrange(-0.14, 0, 7, True):
		S = make_setup("ShortCoil_Cone",r)
		S.ecapCone = [r,r]
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
		
		#BothApertureScan()
		#ECtoWiresScan()
		#ShortCoilScan()
		#MovingCellScan()
		#LengthScan()
		NegApertureScan()
		#DistortionScan()
		#SuperAScan()
		#ConeScan()
		#PosBlockScan()
		
		#GriddingScan()
		
		exit(0)

	outdir = os.environ["ROTSHIELD_OUT"]

	VPP = None

	if options.plot:
	
		#VPP = VarParamPlotter(outdir+"/ShortCoil_MovingCell")
		#VPP.keypos = "br"
		#VPP.setupGraph("Cell center z [m]")
		
		#VPP = VarParamPlotter(outdir+"/ShortCoil_Length")
		#VPP.setupGraph("Coil Length [m]")
		
		#VPP = VarParamPlotter(outdir+"/ShortCoil_Distortion")
		#VPP.keypos = "tr"
		#VPP.setupGraph("Distortion parameter $a$")

		#VPP = VarParamPlotter(outdir+"/ShortCoil_a5")
		#VPP.setupGraph("Distortion parameter $a_3$")
			
		#VPP = VarParamPlotter(outdir+"/ShortCoil_ECdz")
		#VPP.keypos = "bl"
		#VPP.setupGraph("Wires to endcap distance [m]")

		#VPP = VarParamPlotter(outdir+"/ShortCoil_Aperture")
		#VPP.keypos = "bl"
		#VPP.setupGraph("Baseplate aperture [m]")
			
		VPP = VarParamPlotter(outdir+"/ShortCoil_NegAperture")
		VPP.keypos = "bl"
		VPP.setupGraph("Baseplate aperture [m]")

		#VPP = VarParamPlotter(outdir+"/ShortCoil_Gridding")
		#VPP.setupGraph("Endcap grid segments")

		#VPP = VarParamPlotter(outdir+"/ShortCoil_Cone")
		#VPP.setupGraph("Endcap cone offset [m]")

		#VPP = VarParamPlotter(outdir+"/ShortCoil_PosBlock")
		#VPP.setupGraph("Inner disc z offset [m]")
		
		if VPP is not None:
			VPP.makePlot()
			VPP.outputPlot()
			

	if options.fplt:
		if VPP is not None:
			for r,FI in VPP.datlist:
				FI.plotFields(2,0,1,2)
		else:
			
			FI = FieldInfo(outdir+"ShortCoil_PosBlock/X_0.500000/")
			
			FI.plotFields(2,0,1,2) # Bz along z
			FI.plotFields(0,0,2,1) # Bx along y
			FI.plotFields(0,0,1,2) # Bx along z
			
			FI.plotCellProjection(1,0)
			FI.plotCellProjection(1,2)
			FI.plotCellProjection(0,2)


