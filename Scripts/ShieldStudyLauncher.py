#!/usr/bin/python

import os
import time
from optparse import OptionParser
from random import shuffle

class ShieldSection:
	"""Specification for a section of a shield"""
	def __init__(self, mu, cseg, vseg, start="nrd", end="prd"):
		self.mu = mu
		self.cseg = cseg
		self.vseg = vseg
		self.start = start
		self.end = end
		self.off = [[0,0], [0,0]]

	def make_cmd(self):
		"""command line argument for this section"""
		cmd = "add %f %i %i %s %s"%(self.mu,self.cseg,self.vseg,self.start,self.end)
		offset = (self.off[0][0],self.off[0][1],self.off[1][0],self.off[1][1])
		if offset != (0,0,0,0):
			cmd += " mod %f %f %f %f"%offset
		return cmd

class StudySetup:
	def __init__(self,nm,r):
	
		self.r = r						# parameter being varied
	
		self.N = 15						# number of loops in half of cos theta coil
		self.crad = 0.61				# coil radius
		self.clen = 2.5					# coil length
		self.dist = []					# coil distortion parameters
		self.cend = [None,None]			# coild endcap type
		
		self.shieldR = 0				# shield radius (nominal 0.68; 0 = no shield)
		self.shieldL = 0				# shield length (nominal 2.7; 0 = no shield)
		self.nPhi = 64					# shield phi sections
		self.shsects = []				# list of shield sections
			
		self.measCell = [(0,-.2,-.2),(.2,.2,.2)]	# measurement range
		self.measGrid = (7,11,11)					# measurement grid

		self.name = nm					# output directory name
		self.outfl = (nm+"/X_%f")%r		# individual output file name
		self.logdir = "%s/%s/Logs"%(os.environ["ROTSHIELD_OUT"],nm)
		os.system("mkdir -p "+self.logdir)
		
	def set_csgeom(self,clen,crad,dlen=None,drad=None):
		"""Set coil and shield geometry, given coil length/radius and deltas for the shield"""
		self.clen = clen
		self.crad = crad
		if dlen is not None and drad is not None:
			self.shieldL = clen+2*dlen
			self.shieldR = crad+drad
		else:
			self.shieldL = self.shieldR = 0

	def make_cmd(self):
		runcmd = "meas svgrd 0 run %s"%self.outfl
		
		coilcmd = "coil geom %i %f %f"%(self.N,self.clen,self.crad)
		for (n,x) in enumerate(self.dist):
			coilcmd += " dist %i %f"%(n+1,x)
		for (n,s) in enumerate(["neg","pos"]):
			if self.cend[n]:
				coilcmd += " end %s %s"%(self.cend[n],s)
			
		shieldcmd = "shield"
		if self.shieldR and self.shieldL:
			shieldcmd += " geom %f %f %i"%(self.shieldL,self.shieldR,self.nPhi)
		for s in self.shsects:
			shieldcmd += " "+s.make_cmd()
					
		cellcmd = "cell range %f %f %f"%self.measCell[0]
		cellcmd += " %f %f %f"%self.measCell[1]
		cellcmd += " grid %i %i %i"%self.measGrid
				
		return "cd ..; ./RotationShield %s x  %s x  %s x  %s x x > %s/L%f.txt 2>&1\n"%(cellcmd, coilcmd, shieldcmd, runcmd, self.logdir, self.r)

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


	