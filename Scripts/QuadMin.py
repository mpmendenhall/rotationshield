#!/usr/bin/python

from QFile import *
from LinFitter import *
from polynomial import *
import os


###########
# Sequence:
#	- initial guess center, primary axes (orthogonal matrix S), transformed diagonal coefficients lambda_i, and range parameter epsilon
#	- sample points and perform fit in transformed space
#		-> parallel jobs for each sample point location (offload to helper function)
#		-> center offset + composing orthogonal matrix, new lambdas
#		-> if in good agreement, reduce epsilon and repeat


class FitFootprint(QFile):
	"""Optimized pattern for fit point locations"""
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.pts = [array(p.getFirstV("fit_pt")) for p in self.dat.get("fit_pt",[])]



class PointResult(QFile):
	def __init__(self,fname):
		QFile.__init__(self,fname)
		self.x = [float(x) for x in getItem("input","x").split(",")]
		self.y = getItemF("output","y")
		self.nvars = len(self.x)

class QuadMin:
	"""Multivariate `slow' minimum search"""
	
	def __init__(self,basedir):
		self.basedir = basedir
		os.system("mkdir -p "+basedir)
	
	def analyze_step(self,n):
		"""Analyze results of completed step calculation"""
		
		steppath = self.basepath+"/Step_%i/"%n
		
		# initial guess from config file
		qi = QFile(steppath+"/Step.txt")
		self.x0 = [float(x) for x in qi.getItem("initial","x0").split(",")]
		self.nvars = len(self.x0)
		
		# load collected data points
		datpts = [PointResult(steppath+f+"/Point.txt") for f in os.listdir(steppath) if f[:2]=='P_']
		
		# set up quadratic fitter
		
		
		
		
	def minimizer_step(self,n):
		pass

if __name__ == "__main__":

	for ndim in [1,2,3,4,5,6]:
		print
		print
		FF = FitFootprint("QuadMinConfig/OptQuadFitPts_%i.txt"%ndim)
		for p in FF.pts:
			print p.dot(p),p


