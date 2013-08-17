#!/sw/bin/python2.7

import os

import pyx
from pyx import *
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from math import *

from EDM_IO import *


def Vol3Avg(P,ll,ur):
	"Volume average of a polynomial"
	return P.average(0,ll[0],ur[0]).average(0,ll[1],ur[1]).average(0,ll[2],ur[2])

def GradV(V):
	"Gradient of polynomial vector"
	return [Vi.derivative(a) for (a,Vi) in enumerate(V)]


def VaryCoilParam(outdir,basename,varname):
	
	datlist = [ (float(f[len(basename):].split("_")[1]), outdir+"/"+f) for f in os.listdir(outdir) if f[:len(basename)+1]==basename+"_"]
	
	g=graph.graphxy(width=24,height=16,
		x=graph.axis.lin(title=varname),
		y=graph.axis.lin(title="Cell Average Gradients [$\\mu$G/cm]",min=-10,max=10),
		key = graph.key.key(pos="bl",columns=2))
	g.texrunner.set(lfs='foils17pt')
	
	# sample cell dimensions, normal and rotated
	cells = [ [(0.05,-0.05,-0.20),(0.125,0.05,0.20)], [(0.05,-0.20,-0.05),(0.125,0.20,0.05)] ]
	longaxis = [2,1]
	axes = ["x","y","z"]
	# target B0, mG
	B0targ = 30.
		
	for (celln,(cell_ll,cell_ur)) in enumerate(cells):
	
		cellcenter = [(cell_ll[a]+cell_ur[a])*0.5 for a in range(3)]
		print
		print "*** Cell:",cell_ll,cell_ur
		
		# collect data
		gdat = []
		for (r,f) in datlist:
			fitf = open(f+"/Fields/Fieldstats.txt","r")
			B = [read_polynomial(fitf) for a in range(3)]
			# Field average
			Bavg = [Vol3Avg(Bi,cell_ll,cell_ur) for Bi in B]
			# Gradient average
			GradBAvg = [Vol3Avg(Bi,cell_ll,cell_ur) for Bi in GradV(B)]
			# average gradient in uG/cm
			GradScaled = [G*0.01*1000*B0targ/Bavg[0] for G in GradBAvg]
			# sqrt(<dBz/dx^2>) in uG/cm
			dBxdz = B[0].derivative(longaxis[celln])
			rms = sqrt(Vol3Avg(dBxdz*dBxdz,cell_ll,cell_ur))*0.01*1000*B0targ/Bavg[0]
			print r,"B0 =",Bavg,"\tBgrad =",GradScaled,"\tRMS =",rms
			if r>0.3:
				gdat.append([r,GradScaled[0],GradScaled[1],GradScaled[2],rms])
		gdat.sort()
	
		axiscols = {"x":rgb.red,"y":rgb.green,"z":rgb.blue,"S":rgb.black}
		invstyle = [[],[style.linestyle.dashed]][celln]
		invlabel = ["cell along axis","cell rotated"][celln]
		for (n,a) in enumerate(axes):
			g.plot(graph.data.points(gdat,x=1,y=2+n,title="$\\langle \\partial B_{%s} / \\partial %s \\rangle$ "%(a,a)+invlabel),
					[graph.style.line(lineattrs=[axiscols[a],style.linewidth.THick]+invstyle),])
		g.plot(graph.data.points(gdat,x=1,y=5,title="${\\langle (\\partial B_x / \\partial %s)^2 \\rangle}^{1/2}$ "%axes[longaxis[celln]]+invlabel),
					[graph.style.line(lineattrs=[style.linewidth.THick]+invstyle),])

	g.writePDFfile(outdir + "/Plot_%s.pdf"%(basename.split("/")[-1]))


if __name__=="__main__":
	outdir = os.environ["ROTSHIELD_OUT"]

	#VaryCoilParam(outdir+"/Bare_VarLength","CLen","Coil Length [m]")
	#VaryCoilParam(outdir+"/Shielded_VarLength","CLen","Coil Length [m]")
	
	#VaryCoilParam(outdir+"/Bare_L2.5_VarRad","CRad","Coil Radius [m]")
	VaryCoilParam(outdir+"/Shielded_L2.5_VarRad","CRad","Coil Radius [m]")
	
	#VaryCoilParam(outdir+"/Shielded_L2.5_R.45_VarGap","SGap","Shield-to-coil gap [m]")

	#VaryCoilParam(outdir+"/Shielded_L2.5_R.40_VarA","CA","Coil distortion `a'")
