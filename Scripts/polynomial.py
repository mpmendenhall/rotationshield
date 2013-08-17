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
	"""Class for manipulating polynomials in N variables"""
	
	def __init__(self,N):
		self.N = N
		self.C0 = (0,)*self.N
		self.coeffs = {}
		self.varnames = ( "x","y","z","t", "w","v", "u", "a", "b", "c" )
		
	def __add__(self,other):
		"""Add polynomials"""
		p = self.copy()
		if isinstance(other,polynomial):
			assert other.N == self.N
			for k in other.coeffs.keys():
				p.coeffs[k] = p.coeffs.get(k,0) + other.coeffs[k]
		else:
			if other:
				if self.coeffs.has_key(self.C0):
					self.coeffs[self.C0] += other
				else:
					self.coeffs[self.C0] = other
		return p
	
	def __mul__(self,other):
		"""Multiply polynomial or scalar type"""
		p = polynomial(self.N)
		if isinstance(other,polynomial):
			assert other.N == self.N
			for k1 in self.coeffs.keys():
				for k2 in other.coeffs.keys():
					k3 = tuple([k1[i]+k2[i] for i in range(len(k1))])
					p.coeffs[k3] = p.coeffs.get(k3,0) + self.coeffs[k1]*other.coeffs[k2]
		else:
			for k1 in self.coeffs.keys():
				p.coeffs[k1] = self.coeffs[k1]*other 
		return p
		
	def __neg__(self):
		"""Negate polynomial"""
		p = polynomial(self.N)
		for k in self.coeffs.keys():
			p.coeffs[k] = -self.coeffs[k]
		return p
		
	def __sub__(self,other):
		"""Subtract polynomials"""
		return self + (-other)
		
	#def canonical_polyquadratic_vars(self,n):
	#	"""Quadratic polynomial in n variables"""
	#	l = [(0,)*n,]
	#	for i in range(n):
	#		t = [0,]*n
	#		t[i] += 1
	#		l += [tuple(t),]
	#	for i in range(n):
	#		for j in range(i+1):
	#			t = [0,]*n
	#			t[i] += 1
	#			t[j] += 1
	#			l += [tuple(t),]
	#	return l
	
	#def build_polyquadratic(self,n,coefflist = None):
	#	l = self.canonical_polyquadratic_vars(n)
	#	self.coeffs = {}
	#	for i in range(len(l)):
	#		if coefflist:
	#			self.coeffs[l[i]] = float(coefflist[i])
	#		else:
	#			self.coeffs[l[i]] = 1.0
	#	return self

	def eval(self,varvals):
		"""Evaluate at given variable values"""
		s = 0
		for k in self.coeffs.keys():
			ss = self.coeffs[k]
			for i in range(len(k)):
				ss *= varvals[i]**k[i]
			s += ss
		return s
	
	def evalTerms(self,varvals):
		"""Evaluate term-by-term; useful for linear fits"""
		return [ product([varvals[n]**v for (n,v) in enumerate(k)])*self.coeffs[k] for k in self.coeffs.keys() ]
	
	def rescale(self,scalefactors):
		"""Scale each variable by given scale factor x->s*x"""
		for k in self.coeffs.keys():
			self.coeffs[k] *= product([ scalefactors[n]**p for (n,p) in enumerate(k) ])
	
	def evalOne(self,n,xval):
		"""Evaluate out one variable at given value; return scalar if is scalar"""
		p = polynomial(self.N-1)
		for k in self.coeffs.keys():
			nkey = k[:n]+k[n+1:]
			p.coeffs[nkey] = p.coeffs.get(nkey,0) + self.coeffs[k]*(xval**k[n])
		if p.N == 0:
			return p.coeffs.get((),0)
		return p
		
	def derivative(self,n):
		"""Derivative over variable n"""
		p = polynomial(self.N)
		for k in self.coeffs.keys():
			if k[n] == 0:
				continue
			c = self.coeffs[k]*k[n]
			kn = list(k); kn[n] -= 1; kn = tuple(kn);
			p.coeffs[kn] = p.coeffs.get(kn,0) + c
		return p
	
	def indefIntegral(self,n):
		"""Indefinite integral in variable n"""
		p = polynomial(self.N)
		for k in self.coeffs.keys():
			c = self.coeffs[k]/(1.0+k[n])
			kn = list(k); kn[n] += 1; kn = tuple(kn);
			p.coeffs[kn] = p.coeffs.get(kn,0) + c
		return p
		
	def integral(self,n,a,b):
		"""Definite integral integrating out variable n"""
		p = self.indefIntegral(n)
		return p.evalOne(n,b) - p.evalOne(n,a)
		
	def average(self,n,a,b):
		"""Average over range in variable b"""
		return self.integral(n,a,b)*(1.0/(b-a))
		
		
	def copy(self):
		"""Copy into another polynomial"""
		p = polynomial(self.N)
		p.coeffs = self.coeffs.copy()
		p.varnames = self.varnames
		return p
		
	def quadratic_extrema(self):
		lhs,rhs = self.quadratic_derivs_matrix()
		return array(linalg.solve(lhs,rhs).transpose())[0]
		
	def quadratic_derivs_matrix(self):
		lhs = matrix([ [ self.derivative(dv).coeffs[basisv(self.N,v)] for v in range(self.N) ] for dv in range(self.N) ])
		rhs = matrix([ self.derivative(dv).coeffs[(0,)*self.N] for dv in range(self.N) ]).transpose()
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
		"""Display as a string"""
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