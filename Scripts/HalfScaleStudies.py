#!/sw/bin/python2.7

from StudyPlotter import *
from ShieldStudyLauncher import *


def make_halfscale_base(nm, r,
						kdist = 0.005,
						gridt = 58):

	S = StudySetup(nm,r)
	S.nPhi = 64
	S.fields.append(CoilSpec(0.324, 2.146, 15))	# r=1, l=2.5, N=15 cos theta coil
	S.fields[-1].kdist = kdist
	
	sr = 0.02
	srad = 0.362 + sr
	shz = 0.5*2.16
	
	S.shields.append(ShieldSpec(10000, gridt, [(-shz, srad),(shz, srad)], sr, 0.80))
	
	S.measGrid = (9,9,13)
	S.measCell = [(-0.10, -0.10, 0), (0.10, 0.10, 1.0)]

	return S

make_setup = make_halfscale_base
stname = "HS_Base"

def kDistScan():
	SS = StudyScan()
	for (n,r) in  enumerate(unifrange(0, .01, 16)):
		S = make_setup(stname+"/kDist", r, kdist = r)
		S.solfl = "../SC"
		if n == -1:
			os.system(S.make_cmd("RotationShield_Vis"))
		else:
			SS.fsimlist.write(S.make_cmd())
	SS.run()

def basic_halfscale():
	S = make_setup(stname+"/w", 0)
	S.solfl = "../SC"
	S.measGrid = (13,13,51)
	S.measCell = [(-0.12, -0.12, 0), (0.12, 0.12, 1.0)]
	os.system(S.make_cmd("RotationShield_Vis"))

if __name__=="__main__":

	basic_halfscale(); exit(0)
	
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
	
		kDistScan()

	outdir = os.environ["ROTSHIELD_OUT"]

	VPP = None

	if options.plot:

		VPP = VarParamPlotter(outdir+"/"+stname+"/kDist")
		VPP.keypos = "tc"
		VPP.setupGraph("coil distortion `$k$'")

		if VPP is not None:
			VPP.makePlot()
			VPP.outputPlot()

	if options.fplt:
		if VPP is not None:
			for r,FI in VPP.datlist:
				FI.plotFields(2,0,1,2)
		else:
			
			FI = FieldInfo(outdir+"/"+stname+"/kDist/X_0.000000/")
			
			if 1:
				FI.BC.plotFields(2,0,1,2) # Bz along z
				FI.BC.plotFields(0,0,2,1) # Bx along y
				FI.BC.plotFields(0,0,1,2) # Bx along z

