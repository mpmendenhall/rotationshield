#!/usr/bin/python

import numpy
from numpy import *
import numpy.linalg as linalg

def product(l):
	p = 1
	for i in l:
		p *= i
	return p
	
def basisv(n,m=-1):
	l = [0,]*n
	if m>0:
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

	def __call__(self,varvals):
		"""Evaluate at given variable values"""
		s = 0
		for k in self.coeffs.keys():
			ss = self.coeffs[k]
			for i in range(len(k)):
				ss *= varvals[i]**k[i]
			s += ss
		return s
	
	
	def getTerms(self,unityCoeffs=True):
		"""Get each term, optionally with coefficients set to unity"""
		ks = self.coeffs.keys()
		ks.sort()
		Ts = []
		for k in ks:
			P = polynomial(self.N)
			P.coeffs[k] = 1 if unityCoeffs else self.coeffs[k]
			Ts.append(P)
		return Ts
	
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
		
	def quadratic_derivs_matrix(self):
		"""Matrix and vector to solve for minima of quadratic, LHS*x0 = RHS. Note, LHS is symmetric, hence diagonalizable with orthogonal basis."""
		lhs = matrix([ [ self.derivative(dv).coeffs[basisv(self.N,v)] for v in range(self.N) ] for dv in range(self.N) ])
		rhs = matrix([ self.derivative(dv).coeffs[self.C0] for dv in range(self.N) ]).transpose()
		return (lhs, rhs)
	
	def quadratic_extrema(self):
		"""Extremum location for multi-quadratic polynomial"""
		lhs,rhs = self.quadratic_derivs_matrix()
		return array(linalg.solve(lhs,rhs).transpose())[0]
		
	def max_error_bounds(self):
		P = self.copy()
		P.coeffs[(0,)*self.nvars()] = 0
		for i in range(self.nvars()):
			P.coeffs[basisv(self.nvars(),i)] = 0
		m = linalg.pinv(P.quadratic_derivs_matrix()[0])
		return [ m[r,r]/sqrt(2*P( [ -x for x in m.transpose()[r].tolist()[0] ] )) for r in range(self.nvars()) ]
			
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


def lowTriangTerms(nVars, order):
	"""Polynomial of sepcified order and number of variables, e.g. (2,2) -> x^2+xy+y^2+x+y+1"""
	assert order>=0
	if order==0:
		P=polynomial(nVars)
		P.coeffs[P.C0] = 0
		return P
	P = lowTriangTerms(nVars,order-1)
	for n in range(nVars):
		Q = polynomial(nVars)
		Q.coeffs[basisv(nVars,n)] = 0
		P = P+P*Q
	return P


