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
	return P.average(ll[0],ur[1]).average(ll[1],ur[1]).average(ll[2],ur[2])

def GradV(V):
	"Gradient of polynomial vector"
	return [Vi.derivative(a) for (a,Vi) in enumerate(V)]

def ShortCoilVaryRadius():
	outdir = os.environ["ROTSHIELD_OUT"]
	datlist = [ (float(f.split("_")[1]), outdir+"/"+f) for f in os.listdir(outdir) if f[:8]=="HalfCoil"]
	
	# sample cell dimensions
	cell_ll = (0.05,-0.20,-0.05)
	cell_ur = (0.125,0.20,0.05)
	cellcenter = [(cell_ll[a]+cell_ur[a])*0.5 for a in range(3)]
	# target B0, mG
	B0targ = 30.
	
	# collect data
	gdat = []
	for (r,f) in datlist:
		fitf = open(f+"/Fields/Fieldstats.txt","r")
		B = [read_polynomial(fitf) for a in range(3)]
		# Field average
		Bavg = [Vol3Avg(Bi,cell_ll,cell_ur) for Bi in B]
		# Gradient average
		GradBAvg = [Vol3Avg(Bi,cell_ll,cell_ur) for Bi in GradV(B)]
		# average gradient in mG/cm
		GradScaled = [G*0.01*B0targ/Bavg[0] for G in GradBAvg]
		print r,Bavg,GradScaled
		gdat.append([r,GradScaled[0],GradScaled[1],GradScaled[2],sum(GradScaled)])
	gdat.sort()
	
	# plot results
	g=graph.graphxy(width=24,height=16,
		x=graph.axis.lin(title="Coil Radius [m]"),
		y=graph.axis.lin(title="Cell Average Gradients [mG/cm]"),
		key = graph.key.key(pos="bl"))
	g.texrunner.set(lfs='foils17pt')
	
	axes = ["x","y","z"]
	axiscols = {"x":rgb.red,"y":rgb.green,"z":rgb.blue,"S":rgb.black}
	for (n,a) in enumerate(axes):
		g.plot(graph.data.points(gdat,x=1,y=2+n,title="$\\langle \\partial B_{%s} / \\partial %s \\rangle$"%(a,a)), [graph.style.line(lineattrs=[axiscols[a],style.linewidth.THick]),] )
			
	g.writePDFfile(outdir + "/Plot_ShortCoil_Radius.pdf")



if __name__=="__main__":
	ShortCoilVaryRadius()

