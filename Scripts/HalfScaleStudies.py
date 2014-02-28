#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *

def make_halfscale_base(nm,r):

	S = StudySetup(nm,r)
	S.set_csgeom(clen = 2.146, crad = 0.324, dlen = 0.0065, drad = 0.038)
	
	# SC bottom plate
	if 1:
		S.shsects.append(ShieldSection(0,6,12,"nax","nrd"))
		S.shsects[0].off[0] = [-0.051, 0] # 0.102]
		S.shsects[0].off[1] = [-0.051, 0.051]
	
	# main shield
	S.shsects.append(ShieldSection(10000,10,20,"nrd","prd"))
	
	S.measGrid = (7,7,13)
	S.measCell = [(-.18,-0.18,-0.4), (0.18,0.18,0.4)]
	S.dist = [-0.0126]
	
	return S


def OpenEnded():
	SS = StudyScan()
	for x in unifrange(0, 0, 1, True):
		S = make_halfscale_base("Halfscale_Open",x)
		S.shsects[0].off[0][1] = x
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def NegApertureScan():
	SS = StudyScan()
	for x in unifrange(0, 0.2, 7, True):
		S = make_halfscale_base("Halfscale_NegAperture",x)
		S.shsects[0].off[0][1] = x
		SS.fsimlist.write(S.make_cmd())
	SS.run()

def EndcapDistance():
	SS = StudyScan()
	for x in unifrange(-0.05, 0, 7, True):
		S = make_halfscale_base("Halfscale_EC_Dist",x)
		S.shsects[0].off[0][0] = S.shsects[0].off[1][0] = x
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
	
		#OpenEnded()
		#NegApertureScan()
		EndcapDistance()
		
		exit(0)

	outdir = os.environ["ROTSHIELD_OUT"]

	VPP = None

	if options.plot:
	
		VPP = VarParamPlotter(outdir+"/Halfscale_NegAperture")
		VPP.keypos = "br"
		VPP.setupGraph("Bottom plate aperture size [m]")

		if VPP is not None:
			VPP.makePlot()
			VPP.outputPlot()
			

	if options.fplt:
		if VPP is not None:
			for r,FI in VPP.datlist:
				FI.plotFields(2,0,1,2)
		else:
		
			pltdirs = ["Halfscale_NegAperture/X_0.000000/", "Halfscale_Open/X_0.000000/", "Halfscale_NegAperture/X_0.100000/", "Halfscale_EC_Dist/X_0.000000/"]
			
			for d in pltdirs:
				FI = FieldInfo(outdir+d)
				
				FI.BC.ll = [-10,-10,-40]
				FI.BC.ur = [10,10,40]
		
				FI.plotFields(2,0,1,2) # Bz along z
				FI.plotFields(0,0,1,2) # Bx along z
				
				#FI.plotFields(0,0,2,1) # Bx along y
				
				#FI.plotCellProjection(1,0)
				#FI.plotCellProjection(1,2)
				#FI.plotCellProjection(0,2)


