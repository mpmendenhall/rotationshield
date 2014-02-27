#!/sw/bin/python2.7

import os

import pyx
from pyx import *
from pyx.color import rgb, hsb
from pyx.graph.style import symbol
from math import *

from EDM_IO import *
from LinFitter import *
from ArrowPlotter import *

from QFile import *

def rainbow(n):
	return [ hsb((1.0*x)/n,1,1) for x in range(n) ]

class cell_geom(KVMap):
	"""Data cell geometry"""
	def __init__(self,m=None):
		if m:
			self.dat = m.dat
			self.ll = self.getFirstV("ll")
			self.ur = self.getFirstV("ur")
			self.loadFloats(["nx","ny","nz"])

class GeomInfo_File(QFile):
	"""Geometry info file"""
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.cell = cell_geom(self.getFirst("cell"))

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
	def avgB(self):
		return [Vol3Avg(self.B[i],self.ll,self.ur) for i in range(3)]
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

def longAxis(cell):
	"""Determine longest axis of a measurement cell"""
	a = [abs(cell[1][i]-cell[0][i]) for i in range(3)]
	return a.index(max(a))

def linear_zerocrossings(pts):
	"""Linearly interpolate positions of zero-crossings in a list of x,y points"""
	pts.sort()
	zcross = []
	for i in range(len(pts)-1):
		x1,y1 = pts[i]
		x2,y2 = pts[i+1]
		if y1*y2 <= 0:
			zcross.append( [(y2*x1-y1*x2)/(y2-y1),(y2-y1)/(x2-x1)] )
	return zcross

def otherAxis(xi,xj):
	assert xi != xj
	for i in range(3):
		if i not in (xi,xj):
			return i
	return None

class FieldInfo:
	"""Output fields information"""
	def __init__(self,basepath):
		self.basepath = basepath
		self.BC = BCell(basepath+"/Fields/Fieldstats.txt")
		self.GIF = GeomInfo_File(basepath+"/Fields/GeomInfo.txt")

		# target B0, mG
		self.B0targ = 30.
		
		self.BC.ll,self.BC.ur = self.GIF.cell.ll,self.GIF.cell.ur
		self.BC.m2cm()
		self.Bavg = self.BC.normB0(self.B0targ)
		
	def gridLines(self,x0,xi,xj,xk,nptsijk=(3,3,50)):
		"""Bx0 along lines points on axes (xi,xj), [ xk, Bx0(x) ]"""
		self.slicepts = []
		dlines = {}
		for (ni,p_i) in enumerate(unifrange(self.BC.ll[xi],self.BC.ur[xi],nptsijk[0])):
			for (nj,p_j) in enumerate(unifrange(self.BC.ll[xj],self.BC.ur[xj],nptsijk[1])):
				spoints = [ (p_i,p_j,p_k) for p_k in unifrange(self.BC.ll[xk],self.BC.ur[xk],nptsijk[2]) ]
				xpts = [ [0,0,0] for p in spoints ]
				for (n,s) in enumerate(spoints):
					xpts[n][xi] = s[0]
					xpts[n][xj] = s[1]
					xpts[n][xk] = s[2]
				p0 = (xpts[0][xi],xpts[0][xj])
				self.slicepts.append( ((ni,nj), p0) )
				dlines[p0] = [(x[xk],self.BC.B[x0](x)) for x in xpts]
		return dlines
		
	def plotFields(self,x0,xi,xj,xk,nptsijk=(3,3,50)):
		
		# Bx0 slice along fixed xi,xj, varying xk="z"
		
		#x0,xi,xj,xk = 0,0,2,1	# B_x along y
		#x0,xi,xj,xk = 0,0,1,2	# B_x along z
		
		istyles = [[rgb.red],[rgb.green],[rgb.blue]]
		jstyles = [[style.linestyle.dotted,style.linewidth.THick],[style.linestyle.solid,style.linewidth.THick],[style.linestyle.dashed,style.linewidth.Thick]]
		
		g=graph.graphxy(width=24,height=16,
			x=graph.axis.lin(title="$%s$ position [cm]"%axes[xk]),
			y=graph.axis.lin(title="$B_{%s}$ [mG]"%axes[x0]),
			key = graph.key.key(pos="bc",columns=3))
		g.texrunner.set(lfs='foils17pt')
		
		dlines = self.gridLines(x0,xi,xj,xk,nptsijk)		
		for (ni,nj),(p_i,p_j) in self.slicepts:
			g.plot(graph.data.points(dlines[(p_i,p_j)],x=1,y=2,title="$%s,%s = %+.2f,%+.2f$"%(axes[xi],axes[xj],p_i,p_j)),
				[graph.style.line(lineattrs=istyles[ni]+jstyles[nj]),])

		g.writePDFfile(self.basepath + "/Field_B%s_%s.pdf"%(axes[x0],axes[xk]))

	def plotCellProjection(self, xi, xj, nptsijk=None, B0=None, B0scale=None):
	
		if nptsijk is None:
			nptsijk = ( int(self.BC.ur[xi]-self.BC.ll[xi]), int(self.BC.ur[xj]-self.BC.ll[xj]), 5 )
	
		# default to residuals plot
		if B0 is None:
			B0 = self.BC.avgB()
		if B0scale is None:
			B0scale = sqrt(sum([x**2 for x in B0]))*0.005
			
		"""Plot projection of cell fields in xi-xj plane"""
		g=graph.graphxy(width=self.BC.ur[xi]-self.BC.ll[xi]+1, height=self.BC.ur[xj]-self.BC.ll[xj]+1,
			x=graph.axis.lin(title="$%s$ position [cm]"%axes[xi], min=self.BC.ll[xi]-0.5, max=self.BC.ur[xi]+0.5),
			y=graph.axis.lin(title="$%s$ position [cm]"%axes[xj], min=self.BC.ll[xj]-0.5, max=self.BC.ur[xj]+0.5))
		g.texrunner.set(lfs='foils17pt')

		AP = ArrowPlotter(g)
		
		xk = otherAxis(xi,xj)
		dati = self.gridLines(xi,xi,xj,xk,nptsijk)
		datj = self.gridLines(xj,xi,xj,xk,nptsijk)
		
		kcol = rainbow(nptsijk[2])
		
		for nk in range(nptsijk[2]):
			gdat = [ (p[1][0], p[1][1], (dati[p[1]][nk][1]-B0[xi])/B0scale, (datj[p[1]][nk][1]-B0[xj])/B0scale) for p in self.slicepts]
			asty = graph.style.arrow(lineattrs=[style.linewidth.Thick,kcol[nk]])
			AP.plot_arrowdat(gdat,asty,offset=True)

		g.writePDFfile(self.basepath + "/Cell_%s-%s.pdf"%(axes[xi],axes[xj]))
		



class VarParamPlotter:
	"""Plotter for field uniformity varying with parameter"""
	def __init__(self,outdir,basename="X"):

		self.yrange = (None,None)	# y axis plot range
		self.axiscols = {"x":rgb.red,"y":rgb.green,"z":rgb.blue,"S":rgb.black}
		self.keypos = "tl"
		self.g = None
		
		self.outdir = outdir
		self.datlist = [ (float(f[len(basename):].split("_")[1]), FieldInfo(outdir+"/"+f)) for f in os.listdir(outdir) if f[:len(basename)+1]==basename+"_"]
	
	def setupGraph(self,varname):
		
		print "Setting up graph for",varname
		
		self.g = graph.graphxy(width=24,height=16,
			x=graph.axis.lin(title=varname),
			y=graph.axis.lin(title="Cell Average Gradients [$\\mu$G/cm]", min=self.yrange[0], max=self.yrange[1]),
			key = graph.key.key(pos=self.keypos,columns=2))
		self.g.texrunner.set(lfs='foils17pt')
		
		self.g.plot(graph.data.function("y(x)=0",title=None),[graph.style.line(lineattrs=[style.linestyle.dotted,style.linewidth.Thick])])
		
	def makePlot(self,cname="",csty=[],cell=None):
			
			assert self.g is not None
	
			# collect data
			gdat = []
			lax = None
			for (r,FI) in self.datlist:
				if cell:
					FI.BC.ll,FI.BC.ur = cell[0],cell[1]
				lax = longAxis((FI.BC.ll,FI.BC.ur))
				GradScaled = [x*1000 for x in FI.BC.GradBAvg()]
				rms = FI.BC.GradBRMS(0,lax)*1000
				print r,"B0 =",FI.Bavg,"\tBgrad =",GradScaled,"\tRMS =",rms
				gdat.append([r,]+GradScaled+[rms,])
			gdat.sort()

			for (n,a) in enumerate(axes):
				axdat = [(p[0],p[1+n]) for p in gdat]
				print "dB%s/d%s zero crossings: %s"%(a,a,cname),linear_zerocrossings(axdat)
				ptsymb = graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[self.axiscols[a]])
				self.g.plot(graph.data.points(axdat,x=1,y=2,title="$\\langle \\partial B_{%s} / \\partial %s \\rangle$ "%(a,a)+cname),
						[graph.style.line(lineattrs=[self.axiscols[a],style.linewidth.THick]+csty),ptsymb])
			self.g.plot(graph.data.points(gdat,x=1,y=5,title="${\\langle (\\partial B_x / \\partial %s)^2 \\rangle}^{1/2}$ "%axes[lax]+cname),
						[graph.style.line(lineattrs=[style.linewidth.THick]+csty),graph.style.symbol(symbol.circle,size=0.2)])

	def outputPlot(self):
		self.g.writePDFfile(self.outdir + "/FieldUniformity.pdf")
		

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
	
	#FieldPlotter(outdir+"Super_A3/X_-0.100000")
	#FieldPlotter(outdir+"Super_Asearch_4/X_0.018000")
	#FieldPlotter(outdir+"EC_VarRad/X_0.050000")
	
	print FieldPlotter(outdir+"halfscale_open_end/")
	print FieldPlotter(outdir+"halfscale_sc_end/",30*2.76465874521/2.76559952922)



