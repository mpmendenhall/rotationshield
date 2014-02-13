#!/sw/bin/python2.7

import os

import pyx
from pyx import *
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from math import *

from EDM_IO import *
from LinFitter import *

def Vol3Avg(P,ll,ur):
	"Volume average of a polynomial"
	return P.average(0,ll[0],ur[0]).average(0,ll[1],ur[1]).average(0,ll[2],ur[2])

def GradV(V):
	"Gradient of polynomial vector"
	return [Vi.derivative(a) for (a,Vi) in enumerate(V)]

class BCell:
	"""B Field in measurement cell from polynomial"""
	def __init__(self,fname):
		fitf = open(fname,"r")
		self.B = [read_polynomial(fitf) for a in range(3)]
		self.ll = self.ur = ()
	def m2cm(self):
		"""Convert position units from m to cm"""
		for Bi in self.B:
			Bi.rescale([.01,.01,.01])
		self.ll = [100*x for x in self.ll]
		self.ur = [100*x for x in self.ur]
	def normB0(self,Btarg):
		"""Normalize average field to desired quantity"""
		B0avg = Vol3Avg(self.B[0],self.ll,self.ur)
		self.B = [Bi * (Btarg/B0avg) for Bi in self.B]
		return B0avg
	def GradBAvg(self):
		"""Volume average gradient"""
	 	return [Vol3Avg(Bi,self.ll,self.ur) for Bi in GradV(self.B)]
	def GradBRMS(self,xi,xj):
		"""Volume RMS sqrt(<(dBxi/dBxj)^2>)"""
		dBxdz = self.B[xi].derivative(xj)
		return sqrt(Vol3Avg(dBxdz*dBxdz,self.ll,self.ur))

axes = ["x","y","z"]

def VaryCoilParam(outdir,basename,varname):
	
	datlist = [ (float(f[len(basename):].split("_")[1]), outdir+"/"+f) for f in os.listdir(outdir) if f[:len(basename)+1]==basename+"_"]
	
	g=graph.graphxy(width=24,height=16,
		x=graph.axis.lin(title=varname),
		y=graph.axis.lin(title="Cell Average Gradients [$\\mu$G/cm]"), #,min=-7.5,max=7.5),
		key = graph.key.key(pos="tr",columns=2))
	g.texrunner.set(lfs='foils17pt')
	
	# sample cell dimensions, normal and rotated
	cells = [ [(0.05,-0.05,-0.20),(0.125,0.05,0.20)], [(0.05,-0.20,-0.05),(0.125,0.20,0.05)] ]
	longaxis = [2,1]
	# target B0, mG
	B0targ = 30.
		
	for (celln,(cell_ll,cell_ur)) in enumerate(cells):
	
		cellcenter = [(cell_ll[a]+cell_ur[a])*0.5 for a in range(3)]
		print
		print "*** Cell:",cell_ll,cell_ur
		
		# collect data
		gdat = []
		for (r,f) in datlist:
			BC = BCell(f+"/Fields/Fieldstats.txt")
			BC.ll,BC.ur = cell_ll,cell_ur
			BC.m2cm()
			Bavg = BC.normB0(30.)
			GradScaled = [x*1000 for x in BC.GradBAvg()]
			rms = BC.GradBRMS(0,longaxis[celln])*1000
			print r,"B0 =",Bavg,"\tBgrad =",GradScaled,"\tRMS =",rms
			gdat.append([r,]+GradScaled+[rms,])
		gdat.sort()
	
		g.plot(graph.data.function("y(x)=0",title=None),[graph.style.line(lineattrs=[style.linestyle.dotted,style.linewidth.Thick])])
		axiscols = {"x":rgb.red,"y":rgb.green,"z":rgb.blue,"S":rgb.black}
		invstyle = [[],[style.linestyle.dashed]][celln]
		invlabel = ["cell along axis","cell rotated"][celln]
		for (n,a) in enumerate(axes):
			g.plot(graph.data.points(gdat,x=1,y=2+n,title="$\\langle \\partial B_{%s} / \\partial %s \\rangle$ "%(a,a)+invlabel),
					[graph.style.line(lineattrs=[axiscols[a],style.linewidth.THick]+invstyle),])
		g.plot(graph.data.points(gdat,x=1,y=5,title="${\\langle (\\partial B_x / \\partial %s)^2 \\rangle}^{1/2}$ "%axes[longaxis[celln]]+invlabel),
					[graph.style.line(lineattrs=[style.linewidth.THick]+invstyle),])

	g.writePDFfile(outdir + "/Plot_%s.pdf"%(basename.split("/")[-1]))


#
#
def FieldPlotter(basedir, b0Targ=30):

	# load data
	BC = BCell(basedir+"/Fields/Fieldstats.txt")
	
	BC.ll,BC.ur = (0.05,-0.20,-0.05),(0.125,0.20,0.05) # nominal measurement cell
	BC.ll,BC.ur = (-0.15,-0.15,-0.4),(0.15,0.15,0.4) # half-scale measurement region
	
	BC.m2cm()
	Bavg = BC.normB0(b0Targ)
	
	# Bx0 slice along fixed xi,xj, varying xk="z"
	
	#x0,xi,xj,xk = 0,0,2,1	# B_x along y
	x0,xi,xj,xk = 0,0,1,2	# B_x along z
	
	#x0,xi,xj,xk = 1,0,1,2
	
	#x0,xi,xj,xk = 1,0,2,1
	
	nptsi,nptsj,nptsk = 3,3,50
	istyles = [[rgb.red],[rgb.green],[rgb.blue]]
	jstyles = [[style.linestyle.dotted,style.linewidth.THick],[style.linestyle.solid,style.linewidth.THick],[style.linestyle.dashed,style.linewidth.Thick]]
	
	g=graph.graphxy(width=24,height=16,
		x=graph.axis.lin(title="$%s$ position [cm]"%axes[xk]),
		y=graph.axis.lin(title="$B_{%s}$ [mG]"%axes[x0],min=29.9,max=30.1),
		key = graph.key.key(pos="bc",columns=3))
	g.texrunner.set(lfs='foils17pt')

	for (ni,p_i) in enumerate(unifrange(BC.ll[xi],BC.ur[xi],nptsi)):
		for (nj,p_j) in enumerate(unifrange(BC.ll[xj],BC.ur[xj],nptsj)):
			spoints = [ (p_i,p_j,p_k) for p_k in unifrange(BC.ll[xk],BC.ur[xk],nptsk) ]
			xpts = [ [0,0,0] for p in spoints ]
			for (n,s) in enumerate(spoints):
				xpts[n][xi] = s[0]
				xpts[n][xj] = s[1]
				xpts[n][xk] = s[2]
			gdat = [(x[xk],BC.B[x0](x)) for x in xpts]
			g.plot(graph.data.points(gdat,x=1,y=2,title="$%s,%s = %+.2f,%+.2f$"%(axes[xi],axes[xj],p_i,p_j)),
			[graph.style.line(lineattrs=istyles[ni]+jstyles[nj]),])

	g.writePDFfile(basedir + "/Field_B%s_%s.pdf"%(axes[x0],axes[xk]))

	return Bavg


if __name__=="__main__":
	outdir = os.environ["ROTSHIELD_OUT"]

	#VaryCoilParam(outdir+"ShCoilEC_L","X","Coil length [m]")
	#VaryCoilParam(outdir+"Bare_2.4m_OptA","X","Distortion parameter `$a$'")
	
	#FieldPlotter(outdir+"Super_A3/X_-0.100000")
	#FieldPlotter(outdir+"Super_Asearch_4/X_0.018000")
	#FieldPlotter(outdir+"EC_VarRad/X_0.050000")
	
	print FieldPlotter(outdir+"halfscale_open_end/")
	print FieldPlotter(outdir+"halfscale_sc_end/",30*2.76465874521/2.76559952922)



