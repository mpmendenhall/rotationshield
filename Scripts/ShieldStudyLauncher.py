#!/usr/bin/python

import os
import time
from optparse import OptionParser

def Short_Coil_Series(r0,r1,n):
	"""Run simulation series for shielded short coil, varying radius"""
	pcmd = "cd ..; ./RotationShield short %f x\n"
	fsimlist = open("shield_simlist.txt","w")
	for r in [ r0+i*(r1-r0)/(n-1) for i in range(n)]:
		fsimlist.write(pcmd%r)
	fsimlist.close()
	os.system("cat shield_simlist.txt")
	os.system("nice -n 5 parallel -P 6 < shield_simlist.txt")
	os.system("rm shield_simlist.txt")

	
if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option("-k", "--kill", dest="kill", action="store_true", default=False, help="kill running replays")
	parser.add_option("-s", "--short", dest="short", action="store_true", default=False, help="short shield radius sims")
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
	
	if options.short:
		Short_Coil_Series(options.xmin,options.xmax,options.nsteps)
		exit(0)

	print "No option specified. Try --help for listing."
	