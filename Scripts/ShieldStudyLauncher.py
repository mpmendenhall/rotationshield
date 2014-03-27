#!/usr/bin/python

import os
import time
from optparse import OptionParser
from random import shuffle

class ShieldSpec:
	"""Specification for a section of a shield"""
	def __init__(self, mu, nseg, ends, rthick, p=0):
		self.mu = mu			# permeability
		self.nseg = nseg		# segments
		self.rthick = rthick	# end radius (half thickness
		self.ends = ends		# endpoints (1 for slab, 2 for tube
		self.p = p				# adaptive fraction
		
	def make_cmd(self):
		"""command line argument for this section"""
		cmd = ["s","t"][len(self.ends)-1]
		for e in self.ends:
			cmd += " %g %g"%tuple(e)
		cmd += " %g %g %i %g"%(self.rthick,self.mu,self.nseg,self.p)
		return cmd

class CoilSpec:
	"""Specification for cos theta coil"""
	def __init__(self, r=0.61, l=2.5, N=15):
		self.N = N			# half-coil n loops
		self.r = r			# radius
		self.l = l			# length
		self.ends = None	# end wire style; None for default ["arc","arc"]
		self.dist = []		# distortion parameters
	def make_cmd(self):
		cmd = "c geom %i %g %g"%(self.N, self.r, self.l)
		if self.ends:
			cmd += " ends %s %s"%tuple(self.ends)
		for (n,a) in enumerate(self.dist):
			cmd += " dist %i %g"%(n+1,a)
		cmd += " x"
		return cmd

class StudySetup:
	def __init__(self,nm,r):
	
		self.r = r			# parameter being varied
		self.nPhi = 32		# shield phi sections
		self.fields = []	# field sources
		self.shields = []	# material boundaries
		
		self.measCell = [(0,-.2,-.2),(.2,.2,.2)]	# measurement range
		self.measGrid = (7,11,11)					# measurement grid

		self.sng_ep = None

		self.name = nm					# output directory name
		self.outfl = (nm+"/X_%f")%r		# individual output file name
		self.solfl = "none"
		self.logdir = "%s/%s/Logs"%(os.environ["ROTSHIELD_OUT"],nm)
		os.system("mkdir -p "+self.logdir)
		
	def make_cmd(self, rshield="RotationShield"):
		cmd = "dir %s field"%(self.outfl)
		for f in self.fields:
			cmd += " " + f.make_cmd()
		
		cmd += " x bound n %i"%(self.nPhi)
		for s in self.shields:
			cmd += " " + s.make_cmd()
			
		cmd += " x cell range"
		for v in self.measCell:
			for x in v:
				cmd += " %g"%x
		cmd += " grid"
		for x in self.measGrid:
			cmd += " %i"%x
		cmd += " x"
		if self.sng_ep is not None:
			cmd += " ep %g"%self.sng_ep
		cmd += " solve %s meas x"%(self.solfl)
				
		return "cd ..; ./%s %s > %s/L%f.txt 2>&1\n"%(rshield, cmd, self.logdir, self.r)


#####
#####
#####

class StudyScan:
	def __init__(self):
		self.fsimlist = open("shield_simlist.txt","w")		
	def run(self):
		self.fsimlist.close()
		os.system("cat shield_simlist.txt")
		os.system("nice -n 5 parallel < shield_simlist.txt")
		os.system("rm shield_simlist.txt")

#####
#####
#####

def unifrange(xmin,xmax,npts,shuff=False):
	"""Uniformly spaced points, possibly shuffled to random order"""
	if npts == 1:
		return [0.5*(xmin+xmax),]
	l = [ xmin + float(i)/float(npts-1)*(xmax-xmin) for i in range(npts)]
	if shuff:
		shuffle(l)
	return l


	