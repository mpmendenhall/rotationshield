#!/usr/bin/python

import os
import time
from optparse import OptionParser

def Bare_Coil_OptLength(r0,r1,n):
	"""Run simulation series for bare coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 6 11 11 x coil geom 15 %f 0.61 x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,"Bare_VarLength/CLen_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def Shielded_Coil_OptLength(r0,r1,n):
	"""Run simulation series for shielded coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.20 0.20 0.20 grid 6 11 11 x coil geom 15 %f 0.61 x shield geom %f 0.68 x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,r+0.4,"Shielded_VarLength/CLen_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")	

def Bare_Coil_OptRad(r0,r1,n):
	"""Run simulation series for bare coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.175 0.20 0.20 grid 6 11 11 x coil geom 15 2.5 x %f meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,"Bare_L2.5_VarRad/CRad_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

def Shielded_Coil_OptRad(r0,r1,n):
	"""Run simulation series for shielded coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.175 0.20 0.20 grid 6 11 11 x coil geom 15 2.5 %f x shield geom 2.9 %f x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,r+0.07,"Shielded_L2.5_VarRad/CRad_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")	

def Shielded_Coil_OptA(r0,r1,n):
	"""Run simulation series for shielded coil, varying distortion parameter 'a'"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.175 0.20 0.20 grid 6 11 11 x coil geom 15 2.5 0.40 dist 1 %f x shield geom 2.9 0.47 x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(r,"Shielded_L2.5_R.40_VarA/CA_%f"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")	

def Shielded_Coil_OptGap(r0,r1,n):
	"""Run simulation series for shielded coil, varying length"""
	pcmd = "cd ..; ./RotationShield cell range 0 -0.20 -0.20 0.175 0.20 0.20 grid 6 11 11 x coil geom 15 2.5 0.45 x shield geom 2.9 %f x meas svgrd 0 run %s x x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%(0.45+r,"Shielded_L2.5_R.45_VarGap/SGap_%f_m"%r))
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")	

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("--shlen", dest="shlen", action="store_true", default=False, help="variable length shielded coil")
	parser.add_option("--brlen", dest="brlen", action="store_true", default=False, help="variable length bare coil")
	parser.add_option("--shrad", dest="shrad", action="store_true", default=False, help="variable radius shielded coil")
	parser.add_option("--brrad", dest="brrad", action="store_true", default=False, help="variable radius bare coil")
	parser.add_option("--sha", dest="sha", action="store_true", default=False, help="variable distortion shielded coil")
	parser.add_option("--gap", dest="gap", action="store_true", default=False, help="variable coil/shield gap")
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
		Bare_Coil_OptLength(options.xmin,options.xmax,options.nsteps)
		exit(0)

	if options.shrad:
		Shielded_Coil_OptRad(options.xmin,options.xmax,options.nsteps)
		exit(0)

	if options.brrad:
		Bare_Coil_OptRad(options.xmin,options.xmax,options.nsteps)
		exit(0)

	if options.sha:
		Shielded_Coil_OptA(options.xmin,options.xmax,options.nsteps)
		exit(0)

	if options.gap:
		Shielded_Coil_OptGap(options.xmin,options.xmax,options.nsteps)
		exit(0)

	print "No option specified. Try --help for listing."
	