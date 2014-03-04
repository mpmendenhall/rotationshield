#!/sw/bin/python

from StudyPlotter import *
from polynomial import *
from cmath import *
import cmath
from scipy import special

########################
# nEDM geometric phase effects from inhomogeneous fields
# formulae provided by Christopher Swank, Feb. 2014


###################
# Gaussian-cgs units

# physical constants
k_B = 1.3806488e-16 	# Boltzmann contant [cm^2 g / s^2 K]
c =	 29979245800.		# speed of light, [cm / s]

# conversion from SI V/m to Gaussian statV/cm
E_SI_to_GCGS = 1e6/c

###############################
# particle gas properties

class PtclGas:
	"""Definition of particle gas properties"""
	def __init__(self):
		pass

	def p(self, q, w):
		"""Probability density function"""
		# q is a spatial frequency [1/cm]
		# w is an angular velocity [1/s]
		
		# returns units of [s]
		assert False


class PG_3He(PtclGas):
	"""3He superfluid"""
	def __init__(self,T=0.450):
		PtclGas.__init__(self)
		self.gmr = 2.03789e4			# 3He gyromagnetic ratio [1 / Gauss s]
		self.T = T						# gas temperature [K]

	def p(self, q, w):
		m_He3 = 5.00823e-24				# mass of 3He, [g]
		m = 2.4*m_He3					# reduced mass of He3 in superfluid, [g]
		lm = k_B * self.T**8 / (1.6*m)	# collision rate [Hz]

		q = abs(q)
		z = sqrt(m/(2*k_B*self.T))/q
		r = 1j*z*(lm+1j*w)
		x = sqrt(pi)*z*special.wofz(r)
		return 2*x/(1-lm*x)

class PG_n(PtclGas):
	"""UCN gas"""
	def __init__(self,nu_max=350):
		PtclGas.__init__(self)
		self.nu_max = nu_max		# maximum UCN density, [cm/s]
		self.gmr = 1.83247179e4		# neutron gyromagnetic ratio [1 / Gauss s]

	def p(self, q, w):
		q = abs(q)
		x = w/(self.nu_max*q)	# dimensionless constant
		#return 3j*x**2/(2*w)*( (1j*pi-2*atanh(x))*(1/x-x) - 2 )
		return 3j/(2*(self.nu_max*q))*( (1j*pi-2*atanh(x))*(1-x**2) - 2*x) # nicer form for w=0

#############################
# geometric phase calculation

def SBi(BC,i,w,PG):
	"""S_{B_{i}i} [cm^5]*p"""
	
	P = map_poly_to_unit_cell(BC.B[i], BC.ll, BC.ur)
	L = BC.ur[i]-BC.ll[i]	# cell length along axis i, [cm]
	
	s = 0
	for l in range(20)[1::2]:
		FT = Polynomial_Fourier_coeff(P,[(i==a)*0.5*l for a in range(3)])
		pd = PG.p(pi*l/L,w)
		s += (-1)**((l-1)/2) * 4 * L / (pi**2 * l**2) * pd * FT.imag

	return s


def omega_delta_omega(BC, PG, E = 7500000*E_SI_to_GCGS):
	"""Geometric phase frequency shift proportional to E [V/cm]"""
		
	B0 = BC.avgB()[0]	# Average holding field, [Gauss]
	w0 = B0*PG.gmr		# Larmor angular frequency in holding field [Hz]
	
	# subtract off average holding field
	BC.B[0] -= B0
	
	yBy_zBz = BC.B[1]*monomial((0,1,0)) + BC.B[2]*monomial((0,0,1))
	avg_iBi = Vol3Avg(yBy_zBz, BC.ll, BC.ur)	# < x B_x + y B_y > [Gauss cm]
	
	SBxx_SByy = SBi(BC,1,w0,PG) + SBi(BC,2,w0,PG)	# [Gauss cm s]
	
	gm2Ec = PG.gmr**2*E/c	# [statV / Gauss^2 s cm^2] = [1 / Gauss cm s]
	
	dw = (w0*SBxx_SByy.imag/2 - avg_iBi)*gm2Ec	# delta omega [Hz]
	
	return w0,dw

##################
# T_2 calculation

def SBiBi(BC,i,w,PG):
	"""S_{B_i B_i} [Gauss^2 s]"""
	
	P = map_poly_to_unit_cell(BC.B[i], BC.ll, BC.ur)
	Lx,Ly,Lz = [1.0*abs(BC.ur[a]-BC.ll[a]) for a in range(3)]
	
	s = 0
	lmax = 5 # number of terms to calculate
	
	#for lx in range(-lmax,lmax+1):
	#	for ly in range(-lmax,lmax+1):
	#		for lz in range(-lmax,lmax+1):
			
	for lx in range(lmax+1):
		for ly in range(lmax+1):
			for lz in range(lmax+1):

				nzero = sum([l==0 for l in (lx,ly,lz)])	# number of l_i which are 0
				if nzero == 3:
					continue	# exclude B0^2 term
					
				q = sqrt( (lx/Lx)**2 + (ly/Ly)**2 + (lz/Lz)**2 )*pi
				p = PG.p(q,w)
				FT = Polynomial_Fourier_coeff(P,(0.5*lx, 0.5*ly, 0.5*lz))
				FTm = Mirrored_Polynomial_Fourier_coeff(P,(lx, ly, lz))
				
				#s += p*FT*FTm
				s += p * (FT*FTm).real * 2**(3-nzero)
				
	print "S_B%iB%i(%f) ="%(i,i,w),s
		
	return s

def omega_T2i(BC,PG):
	"""1/T_2 dephasing time in inhomogeneous field [Hz]"""
	B0 = BC.avgB()[0]	# Average holding field, [Gauss]
	w0 = B0*PG.gmr		# Larmor angular frequency in holding field [Hz]
	return w0, PG.gmr**2/4. * (2*SBiBi(BC,0,0,PG) + SBiBi(BC,1,w0,PG) + SBiBi(BC,2,w0,PG)).real	# 1/T_2 [Hz]


###########################
# test example calculations

if __name__ == "__main__":
	
	# define magnetic gradients cell
	B0 = [0.030, 0, 0]		# B_0 field [Gauss]
	#gradB = [0, 0, -1e-7] 	# linear field gradient [G/cm]
	gradB = [-1e-7, 0, 0] 	# linear field gradient [G/cm]
	BC = BCell()
	for a in range(3):
		BC.B[a] += B0[a]
		BC.B[a] += monomial(basisv(3,a),gradB[a])
	BC.ll,BC.ur = (-3.75,-20,-5.1),(3.75,20,5.1)	# field cell dimansions [cm]

	print "Calculating phase shift..."
	gdat_n = []
	gdat_3He = []
	for b0 in unifrange(0.0001,0.050,20):
		BC.B[0].coeffs[BC.B[0].C0] = b0
		w0,T2i = omega_T2i(BC,PG_n())
		w0,dw = omega_delta_omega(BC,PG_n())
		gdat_n.append([w0/(2*pi),dw,T2i])
		print "n:",w0,dw,T2i
		
		BC.B[0].coeffs[BC.B[0].C0] = b0
		w0,T2i = omega_T2i(BC,PG_3He())
		w0,dw = omega_delta_omega(BC,PG_3He())
		gdat_3He.append([w0/(2*pi),dw,T2i])
		print "3He:",w0,dw,T2i


	print "Plotting result..."

	if 0:
		gdw = graph.graphxy(width=16,height=12,
			x=graph.axis.lin(title="Larmor frequency $f$ [Hz]",min=80,max=120),
			y=graph.axis.lin(title="Geometric phase $\\delta\\omega$",min=-3e-8,max=3e-8),
			key = graph.key.key())
		gdw.texrunner.set(lfs='foils17pt')

		gdw.plot(graph.data.function("y(x)=0",title=None),[graph.style.line([style.linestyle.dotted]),])
		gdw.plot(graph.data.points(gdat_n,x=1,y=2,title="neutron"),[graph.style.line(),])
		gdw.plot(graph.data.points(gdat_3He,x=1,y=2,title="$^3$He"),[graph.style.line([style.linewidth.THick]),])
		gdw.writePDFfile("/Users/michael/Desktop/GeomPhase.pdf")

	if 1:
		gT2 = graph.graphxy(width=16,height=12,
			x=graph.axis.lin(title="Larmor frequency $f$ [Hz]"),
			y=graph.axis.lin(title="Dephasing $1/T_2$ [Hz]"),
			key = graph.key.key())
		gT2.texrunner.set(lfs='foils17pt')
		gT2.plot(graph.data.points(gdat_n,x=1,y=3,title="neutron"),[graph.style.line(),])
		gT2.plot(graph.data.points(gdat_3He,x=1,y=3,title="$^3$He"),[graph.style.line([style.linewidth.THick]),])
		gT2.writePDFfile("/Users/michael/Desktop/T2.pdf")

