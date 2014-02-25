#!/usr/bin/python

#./RotationShieldVis cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 7 11 11 x coil geom 15 3.0 0.61 x shield geom 3.06 0.68 ecap 15 0.25 0 x meas svgrd 0 run foo x

import os
import time
from optparse import OptionParser
from random import shuffle

class StudySetup:
	def __init__(self,nm,r):
	
		self.r = r						# parameter being varied
	
		self.N = 15						# number of loops in half of cos theta coil
		self.crad = 0.61				# coil radius
		self.clen = 2.5					# coil length
		self.dist = []					# coil distortion parameters
		
		self.shieldR = 0				# shield radius (nominal 0.68; 0 = no shield)
		self.shieldL = 0				# shield length (nominal 2.7; 0 = no shield)
		self.sGrid = (10,20,128)		# shield gridding fixed,variable,phi (10,20,128)
		
		
		self.ecapS = [0,0]				# number of endcap segments (0 = no endcaps)
		self.ecapR = [0,0]				# endcap hole radius
		self.ecapMu = [0,0]				# endcap mu (0 = superconductor)
		self.ecapCone = [0,0]			# endcap cone offset
		
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
		if dlen and drad:
			self.shieldL = clen+2*dlen
			self.shieldR = crad+drad
		else:
			self.shieldL = self.shieldR = 0

	def set_ecaps(self,nSeg,r0,mu,end=None):
		if end is None:
			self.set_ecaps(nSeg,r0,mu,0)
			self.set_ecaps(nSeg,r0,mu,1)
		else:
			self.ecapS[end] = nSeg
			self.ecapR[end] = r0
			self.ecapMu[end] = mu

	def make_cmd(self):
		runcmd = "meas svgrd 0 run %s"%self.outfl
		
		coilcmd = "coil geom %i %f %f"%(self.N,self.clen,self.crad)
		for (n,x) in enumerate(self.dist):
			coilcmd += " dist %i %f"%(n+1,x)
		
		shieldcmd = "shield"
		if self.shieldR and self.shieldL:
			shieldcmd += " geom %f %f"%(self.shieldL,self.shieldR)
			shieldcmd += " grid %i %i %i"%self.sGrid
		for end in [0,1]:
			ename = ["neg","pos"][end]
			if self.ecapS[end] and self.ecapR[end] != self.shieldR:
				shieldcmd += " ecap %i %f %f %s"%(self.ecapS[end],self.ecapR[end],self.ecapMu[end],ename)
			if self.ecapCone[end]:
				shieldcmd += " cone %f %s"%(self.ecapCone[end],ename)
					
		cellcmd = "cell range %f %f %f"%self.measCell[0]
		cellcmd += " %f %f %f"%self.measCell[1]
		cellcmd += " grid %i %i %i"%self.measGrid
				
		return "cd ..; ./RotationShield %s x  %s x  %s x  %s x x > %s/L%f.txt 2>&1\n"%(cellcmd, coilcmd, shieldcmd, runcmd, self.logdir, self.r)


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

def Bare_Misc(r0,r1,n):
	fsimlist = open("shield_simlist.txt","w")
	for r in unifrange(r0,r1,n,True):
		S = StudySetup("Bare_A-.0021_VarL",r)
		S.set_csgeom(r,0.61)
		S.dist = [-0.0021]
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def Shield_Misc(r0,r1,n):
	fsimlist = open("shield_simlist.txt","w")
	for r in unifrange(r0,r1,n,True):
		S = StudySetup("ShCoilEC_L",r)
		#S.set_csgeom(2.5,0.61,.007,0.07)
		#S.dist = [-.0075,0,-.075,0,r]
		S.set_csgeom(r,0.61,0.007,0.07)
		S.dist = [-0.00475]
		S.ecapS = 15
		S.ecapR = 0.10
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")


#####
#####
#####

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("--misc", dest="misc", action="store_true", default=False, help="perform misc study specified in script")
	parser.add_option("--bare", dest="bare", action="store_true", default=False, help="perform bare coil study")
	parser.add_option("--xmin", type="float", dest="xmin", default=0.5)
	parser.add_option("--xmax", type="float", dest="xmax", default=1.5)
	parser.add_option("--nsteps", type="int", dest="nsteps", default=12, help="Number of steps in scan")
	
	options, args = parser.parse_args()
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 RotationShield")
		os.system("killall -9 ShieldStudyLauncher.py")
		exit(0)
	
	if options.misc:
		Shield_Misc(options.xmin,options.xmax,options.nsteps)
		exit(0)
	
	if options.bare:
		Bare_Misc(options.xmin,options.xmax,options.nsteps)
		exit(0)
	
	print "No option specified. Try --help for listing."
	