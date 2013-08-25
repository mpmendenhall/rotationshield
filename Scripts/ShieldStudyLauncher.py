#!/usr/bin/python

#./RotationShieldVis cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 7 11 11 x coil geom 15 3.0 0.61 x shield geom 3.06 0.68 ecap 15 0.25 0 x meas svgrd 0 run foo x

import os
import time
from optparse import OptionParser

class StudySetup:
	def __init__(self,nm,r):
		self.N = 15
		self.r = r
		
		self.crad = 0.61
		self.shieldR = 0.68
		self.clen = 2.5
		self.shieldL = 2.7
		
		self.cellcmd = "cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 7 11 11"
		self.name = nm
		self.outfl = (nm+"/X_%f")%r
		self.dist = []
		self.shieldR = self.shieldL = 0
		self.ecapS = 0
		self.ecapR = 0
		self.ecapMu = 0
		
		
		self.logdir = "%s/%s/Logs"%(os.environ["ROTSHIELD_OUT"],nm)
		os.system("mkdir -p "+self.logdir)
		
	def set_csgeom(self,clen,crad,dlen,drad):
		"""Set coil and shield geometry, given coil length/radius and deltas for the shield"""
		self.clen = clen
		self.shieldL = clen+2*dlen
		self.crad = crad
		self.shieldR = crad+drad
		
	def make_cmd(self):
		runcmd = "meas svgrd 0 run %s"%self.outfl
		coilcmd = "coil geom %i %f %f"%(self.N,self.clen,self.crad)
		for (n,x) in enumerate(self.dist):
			coilcmd += " dist %i %f"%(n+1,x)
		shieldcmd = "shield"
		if self.shieldR and self.shieldL:
			shieldcmd += " geom %f %f"%(self.shieldL,self.shieldR)
		if self.ecapS:
			shieldcmd += " ecap %i %f %f"%(self.ecapS,self.ecapR,self.ecapMu)
			
		return "cd ..; ./RotationShield %s x  %s x  %s x  %s x x >> %s/L%f.txt\n"%(self.cellcmd,coilcmd,shieldcmd,runcmd,self.logdir,self.r)


def ShCoil_ECDist(r0,r1,n):
	"""Shielded coil endcap distance"""
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		S = StudySetup("ShCoil_ECDist",r)
		S.set_csgeom(2.5,0.61,2*r,0.07)
		S.dist = [-0.00475]
		S.ecapS = 15
		S.ecapR = 0.10
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")


def ShCoil_ECRad(r0,r1,n):
	"""Shielded coil endcap distance"""
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		S = StudySetup("ShCoil_ECRad",r)
		S.set_csgeom(2.5,0.61,.007,0.07)
		S.dist = [-0.00475]
		S.ecapS = 15
		S.ecapR = r
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def ShCoilEC_A(r0,r1,n):
	"""Shielded coil endcap distance"""
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		S = StudySetup("ShCoilEC_A",r)
		S.set_csgeom(2.5,0.61,.007,0.07)
		S.dist = [r,]
		S.ecapS = 15
		S.ecapR = 0.10
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def SuperShield_Len(r0,r1,n):
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		S = StudySetup("Super_Len",r)
		S.set_csgeom(r,0.61,.007,0.07)
		S.dist = [-0.00475,]
		S.ecapS = 15
		S.ecapR = 0.10
		fsimlist.write(S.make_cmd())
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")




#
#
#

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("--misc", dest="misc", action="store_true", default=False, help="perform misc study specified in script")
	parser.add_option("--xmin", type="float", dest="xmin", default=0.5)
	parser.add_option("--xmax", type="float", dest="xmax", default=1.5)
	parser.add_option("--nsteps", type="int", dest="nsteps", default=12, help="Number of steps in scan")
	
	options, args = parser.parse_args()
	if options.kill:
		os.system("killall -9 parallel")
		os.system("killall -9 RotationShield")
		os.system("killall -9 ShieldStudyLauncher.py")
		exit(0)
	
	#if len(os.popen("ps -a | grep ShieldStudyLauncher").readlines()) > 1:
	#	print "Already running! I die!"
	#	exit(1)
	
	if options.misc:
		SuperShield_Len(options.xmin,options.xmax,options.nsteps)
		exit(0)
	
	print "No option specified. Try --help for listing."
	