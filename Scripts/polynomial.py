#!/usr/bin/python

import numpy
from numpy import *
import numpy.linalg as linalg

def product(l):
	p = 1
	for i in l:
		p *= i
	return p
	
def basisv(n,m):
	l = [0,]*n
	l[m] = 1
	return tuple(l)


# polynomial, represented by a dictionary with tuples of exponents as the keys
# mapped to the coefficient of the term
class polynomial:

	coeffs = {}
	varnames = ( "x","y","z","t", "w","v", "u", "a", "b", "c" )
	
	def __init__(self):
		self.coeffs = {}
		
	def __add__(self,other):
		p = self.copy()
		for k in other.coeffs.keys():
			p.coeffs[k] = p.coeffs.get(k,0) + other.coeffs[k]
		return p
	
	def __mul__(self,other):
		p = polynomial()
		if isinstance(other,polynomial):
			for k1 in self.coeffs.keys():
				for k2 in other.coeffs.keys():
					k3 = tuple([k1[i]+k2[i] for i in range(len(k1))])
					p.coeffs[k3] = p.coeffs.get(k3,0) + self.coeffs[k1]*other.coeffs[k2]
		else:
			for k1 in self.coeffs.keys():
				p.coeffs[k1] = self.coeffs[k1]*other 
		return p
		
	def __neg__(self):
		p = polynomial()
		for k in self.coeffs.keys():
			p.coeffs[k] = -self.coeffs[k]
		return p
		
	def __sub__(self,other):
		return self + (-other)
	
	def nvars(self):
		return len(self.coeffs.keys()[0]);
		
	def addConst(self,c):
		self.coeffs[(0,)*self.nvars()] += c
		return self
		
	def canonical_polyquadratic_vars(self,n):
		l = [(0,)*n,]
		for i in range(n):
			t = [0,]*n
			t[i] += 1
			l += [tuple(t),]
		for i in range(n):
			for j in range(i+1):
				t = [0,]*n
				t[i] += 1
				t[j] += 1
				l += [tuple(t),]
		return l
	
	def build_polyquadratic(self,n,coefflist = None):
		l = self.canonical_polyquadratic_vars(n)
		self.coeffs = {}
		for i in range(len(l)):
			if coefflist:
				self.coeffs[l[i]] = float(coefflist[i])
			else:
				self.coeffs[l[i]] = 1.0
		return self

	def eval(self,varvals):
		s = 0
		for k in self.coeffs.keys():
			ss = self.coeffs[k]
			for i in range(len(k)):
				ss *= varvals[i]**k[i]
			s += ss
		return s
	
	def evalTerms(self,varvals):
		return [ product([varvals[n]**v for (n,v) in enumerate(k)])*self.coeffs[k] for k in self.coeffs.keys() ]
	
	def rescale(self,scalefactors):
		for k in self.coeffs.keys():
			self.coeffs[k] *= product([ scalefactors[n]**p for (n,p) in enumerate(k) ])
	
	def evalOne(self,xval):
		p = polynomial()
		for k in self.coeffs.keys():
			p.coeffs[k[1:]] = p.coeffs.get(k[1:],0) + self.coeffs[k]*(xval**k[0])
		if () in p.coeffs.keys():
			return p.coeffs[()]
		return p
		
	def derivative(self,n):
		p = polynomial()
		for k in self.coeffs.keys():
			if k[n] == 0:
				continue
			c = self.coeffs[k]*k[n]
			kn = list(k); kn[n] -= 1; kn = tuple(kn);
			p.coeffs[kn] = p.coeffs.get(kn,0) + c
		return p
	
	def indefIntegral(self,n):
		p = polynomial()
		for k in self.coeffs.keys():
			c = self.coeffs[k]/(1.0+k[n])
			kn = list(k); kn[n] += 1; kn = tuple(kn);
			p.coeffs[kn] = p.coeffs.get(kn,0) + c
		return p
		
	def integral(self,a,b,n=0):
		p = self.indefIntegral(n)
		return p.evalOne(b) - p.evalOne(a)
		
	def average(self,a,b,n=0):
		return self.integral(a,b,n)*(1.0/(b-a))
		
		
	def copy(self):
		p = polynomial()
		p.coeffs = self.coeffs.copy()
		p.varnames = self.varnames
		return p
		
	def quadratic_extrema(self):
		lhs,rhs = self.quadratic_derivs_matrix()
		return array(linalg.solve(lhs,rhs).transpose())[0]
		
	def quadratic_derivs_matrix(self):
		n = self.nvars()
		lhs = matrix([ [ self.derivative(dv).coeffs[basisv(n,v)] for v in range(n) ] for dv in range(n) ])
		rhs = matrix([ self.derivative(dv).coeffs[(0,)*n] for dv in range(n) ]).transpose()
		return (lhs, rhs)
	
	def max_error_bounds(self):
		P = self.copy()
		P.coeffs[(0,)*self.nvars()] = 0
		for i in range(self.nvars()):
			P.coeffs[basisv(self.nvars(),i)] = 0
		m = linalg.pinv(P.quadratic_derivs_matrix()[0])
		return [ m[r,r]/sqrt(2*P.eval( [ -x for x in m.transpose()[r].tolist()[0] ] )) for r in range(self.nvars()) ]
	
	def testpoly(self):
			self.coeffs = { (1,0):1, (0,1):2, (1,1):3 }
			
	def tostring(self):
		s = ""
		ks = self.coeffs.keys()
		ks.sort()
		for k in ks:
			if self.coeffs[k] == 0:
				continue
			s += "\t%+g"%self.coeffs[k];
			for i in range(len(k)):
				if k[i] == 0:
					continue
				s += " %s"%self.varnames[i]
				if k[i] > 1:
					s += "^%i"%k[i]
					
		return s[1:]