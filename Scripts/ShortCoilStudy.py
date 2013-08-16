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
	return P.average(ll[0],ur[0]).average(ll[1],ur[1]).average(ll[2],ur[2])

def GradV(V):
	"Gradient of polynomial vector"
	return [Vi.derivative(a) for (a,Vi) in enumerate(V)]

def ShortCoilVaryRadius(basename):
	outdir = os.environ["ROTSHIELD_OUT"]
	datlist = [ (float(f.split("_")[1]), outdir+"/"+f) for f in os.listdir(outdir) if f[:len(basename)+1]==basename+"_"]
	
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
			
	g.writePDFfile(outdir + "/Plot_%s_Radius.pdf"%basename)


def BareCoilVaryLength(basename):
	outdir = os.environ["ROTSHIELD_OUT"]
	datlist = [ (float(f.split("_")[2]), outdir+"/"+f) for f in os.listdir(outdir) if f[:len(basename)+1]==basename+"_"]
	
	g=graph.graphxy(width=24,height=16,
		x=graph.axis.lin(title="Coil Length [m]"),
		y=graph.axis.lin(title="Cell Average Gradients [mG/cm]"),
		key = graph.key.key(pos="br"))
	g.texrunner.set(lfs='foils17pt')
	
	# sample cell dimensions, normal and rotated
	cells = [ [(0.05,-0.05,-0.20),(0.125,0.05,0.20)], [(0.05,-0.20,-0.05),(0.125,0.20,0.05)] ]
	# target B0, mG
	B0targ = 30.
		
	for (celln,(cell_ll,cell_ur)) in enumerate(cells):
	
		cellcenter = [(cell_ll[a]+cell_ur[a])*0.5 for a in range(3)]
		
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
			if r >= 2:
				gdat.append([r,GradScaled[0],GradScaled[1],GradScaled[2],sum(GradScaled)])
		gdat.sort()
	
		axes = ["x","y","z"]
		axiscols = {"x":rgb.red,"y":rgb.green,"z":rgb.blue,"S":rgb.black}
		invstyle = [[],[style.linestyle.dashed]][celln]
		invlabel = ["cell along axis","cell rotated"][celln]
		for (n,a) in enumerate(axes):
			g.plot(graph.data.points(gdat,x=1,y=2+n,title="$\\langle \\partial B_{%s} / \\partial %s \\rangle$ "%(a,a)+invlabel),
					[graph.style.line(lineattrs=[axiscols[a],style.linewidth.THick]+invstyle),])
				
	g.writePDFfile(outdir + "/Plot_%s_Length.pdf"%basename)


if __name__=="__main__":
	#ShortCoilVaryRadius("ShortCoilBare")
	#ShortCoilVaryRadius("ShortCoil")
	BareCoilVaryLength("BareCoil_15")
