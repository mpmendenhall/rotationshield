#!/usr/bin/python

from math import *
from polynomial import polynomial

# return [ [[x,y,z],[bx,by,bz]], ... ] from xyz-grid fieldmap file
def read_xyzgrid(fname):
	lines = [ l for l in open(fname,"r").readlines() if l[0] != '#' ]
	return [ [s[:3],s[3:]] for s in [ [ float(x) for x in l.split()] for l in lines ] ]

# return [ [[x,y,z],bx,dbx,time,current], ... ] from Justin Chen's field map files
def read_chengrid(fname):
	lines = [ l for l in open(fname,"r").readlines() if len(l)>10 and l[0] != '%' ]
	return [ [ [p/6552.0 for p in s[:3]], ] + s[3:] for s in [ [ float(x) for x in l.split()] for l in lines ] ]

# read the next polynomial in the file
def read_polynomial(fin):
	l = fin.readline()
	while len(l) and l[0] != ' ':
		l = fin.readline()
	
	p = polynomial()
	while len(l) and l[0] == ' ':
		t = [ float(x) for x in l.split() ]
		p.coeffs[ tuple(t[1:]) ] = t[0]
		l = fin.readline()
		
	return p
	
# center point of sampling grid	
def centerpoint(fgrid):
	# all points at y=z=0
	c = [ p for p in fgrid if p[0][1] == 0 and p[0][2] == 0 ]
	# middle point of these
	return c[len(c)/2]
	
# convert gridded field to deviation relative to B0x
def field_to_B0x_deviation(fgrid):
	B0x = centerpoint(fgrid)[1][0]
	return [ [ p[0], [(p[1][0]*1.0 - B0x)/B0x,p[1][1]/B0x,p[1][2]/B0x] ] for p in fgrid ]
	

# split gridded field into a set of lines parallel to the given axis	
def axis_parallel_lines(fgrid,mainaxis):
	otheraxes = [ i for i in range(3) if i != mainaxis ]
	lines = {}
	for p in fgrid:
		k = tuple([ p[0][i] for i in otheraxes ])
		lines[k] = lines.get(k,[]) + [p,]
	return lines