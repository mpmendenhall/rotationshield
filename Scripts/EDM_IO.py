#!/usr/bin/python

from math import *
from polynomial import polynomial

# return [ [[x,y,z],[bx,by,bz]], ... ] from xyz-grid fieldmap file
def read_xyzgrid(fname):
    lines = [ l for l in open(fname,"r").readlines() if l[0] != '#' ]
    return [ [s[:3],s[3:]] for s in [ [ float(x) for x in l.split()] for l in lines ] ]

# read the next polynomial in the file
def read_polynomial(fin):
    l = fin.readline()
    while len(l) and l[0] != ' ':
        l = fin.readline()
    
    p = polynomial(3)
    while len(l) and l[0] == ' ':
        t = [ float(x) for x in l.split() ]
        p.coeffs[ tuple(t[1:]) ] = t[0]
        l = fin.readline()
        
    return p

# split gridded field into a set of lines parallel to the given axis    
def axis_parallel_lines(fgrid,mainaxis):
    otheraxes = [ i for i in range(3) if i != mainaxis ]
    lines = {}
    for p in fgrid:
        k = tuple([ p[0][i] for i in otheraxes ])
        lines[k] = lines.get(k,[]) + [p,]
    return lines