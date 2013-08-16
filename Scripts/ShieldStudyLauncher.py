#!/usr/bin/python

import os
import time
from optparse import OptionParser

def Bare_Coil_Optlength(r0,r1,n):
	"""Run simulation series for bare coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 6 11 11 x coil 15 %f 0.61 meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,"Bare_VarLength/CLen_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def Shielded_Coil_OptLength(r0,r1,n):
	"""Run simulation series for shielded coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 6 11 11 x coil 15 %f 0.61 shield geom %f 0.68 x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,r+0.4,"Shielded_VarLength/CLen_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")	

	
if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("--shlen", dest="shlen", action="store_true", default=False, help="variable length shielded coil")
	parser.add_option("--brlen", dest="brlen", action="store_true", default=False, help="variable length bare coil")
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
	
	if options.shlen:
		Shielded_Coil_OptLength(options.xmin,options.xmax,options.nsteps)
		exit(0)

	if options.brlen:
		Bare_Coil_Optlength(options.xmin,options.xmax,options.nsteps)
		exit(0)

	print "No option specified. Try --help for listing."
	