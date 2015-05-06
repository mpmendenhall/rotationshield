#!/usr/bin/python

from StudyPlotter import *
from polynomial import *
from math import *
import numpy
import cmath
from scipy import special
import os
import time

###################################################
# spin precession effects from inhomogeneous fields
# formalism provided by Christopher Swank, March 2014
#
# everything in Gaussian-cgs units

# physical constants
k_B = 1.3806488e-16     # Boltzmann contant [cm^2 g / s^2 K]
c =     29979245800.        # speed of light, [cm / s]

# conversion from SI E=V/m to Gaussian-cgs statVolts/cm
E_SI_to_GCGS = 1e6/c


###############################
# particle gas properties
# provides gyromagnetic ratio
# and p(q,omega)

class PtclGas:
    """Definition of particle gas properties"""
    def __init__(self):
        self.gmr = None    # gyromagnetic ratio  [1 / Gauss s]

    def p(self, q, w):
        """Probability density function"""
        # q is a spatial frequency [1/cm]
        # w is an angular velocity [Hz]
        # returns units of [s]
        assert False

class PG_3He(PtclGas):
    """3He superfluid"""
    def __init__(self,T=0.450):
        PtclGas.__init__(self)
        self.gmr = 2.03789e4            # 3He gyromagnetic ratio [1 / Gauss s]
        self.T = T                        # gas temperature [K]
        self.name = "$^3$He"
        
    def p(self, q, w):
        m_He3 = 5.00823e-24                # mass of 3He, [g]
        m = 2.4*m_He3                    # reduced mass of He3 in superfluid, [g]
        lm = k_B * self.T**8 / (1.6*m)    # collision rate [Hz]
    
        q = abs(q)
        z = sqrt(m/(2*k_B*self.T))/q
        r = 1j*z*(lm+1j*w)
        x = sqrt(pi)*z*special.wofz(r)
        
        return 2*x/(1-lm*x)

class PG_n(PtclGas):
    """UCN gas"""
    def __init__(self,nu_max=350):
        PtclGas.__init__(self)
        self.nu_max = nu_max        # most probably UCN velocity [cm/s]
        self.gmr = 1.83247179e4        # neutron gyromagnetic ratio [1 / Gauss s]
        self.name = "neutron"
        
    def p(self, q, w):
        q = abs(q)
        x = w/(self.nu_max*q)    # dimensionless constant
        #return 3j*x**2/(2*w)*( (1j*pi-2*atanh(x))*(1/x-x) - 2 )
        return 3j/(2*(self.nu_max*q))*( (1j*pi-2*cmath.atanh(x))*(1-x**2) - 2*x) # nicer form for w=0




#############################
# geometric phase calculation

def SBi(BC,i,w,PG):
    """S_{B_{i}i} [cm^5]*p"""
    
    P = map_poly_to_unit_cell(BC.B[i], BC.ll, BC.ur)
    L = BC.ur[i]-BC.ll[i]    # cell length along axis i, [cm]
    
    s = 0
    for l in range(20)[1::2]:
        FT = Polynomial_Fourier_coeff(P,[(i==a)*0.5*l for a in range(3)])
        pd = PG.p(pi*l/L,w)
        s += (-1)**((l-1)/2) * 4 * L / (pi**2 * l**2) * pd * FT.imag

    return s


def omega_delta_omega(BC, PG, E = 7500000*E_SI_to_GCGS):
    """Geometric phase frequency shift proportional to E [V/cm]"""
        
    B0 = BC.avgB()[0]    # Average holding field, [Gauss]
    w0 = B0*PG.gmr        # Larmor angular frequency in holding field [Hz]
    
    # subtract off average holding field
    BC.B[0] -= B0
    
    yBy_zBz = BC.B[1]*monomial((0,1,0)) + BC.B[2]*monomial((0,0,1))
    avg_iBi = Vol3Avg(yBy_zBz, BC.ll, BC.ur)    # < x B_x + y B_y > [Gauss cm]
    
    SBxx_SByy = SBi(BC,1,w0,PG) + SBi(BC,2,w0,PG)    # [Gauss cm s]
    
    gm2Ec = PG.gmr**2*E/c    # [statV / Gauss^2 s cm^2] = [1 / Gauss cm s]
    
    dw = (w0*SBxx_SByy.imag/2 - avg_iBi)*gm2Ec    # delta omega [Hz]
    
    return w0,dw


    

##################
# T_2 calculation

def SBiBi(BC,i,w,PG):
    """S_{B_i B_i} [Gauss^2 s]"""
    #tstart = time.clock()
    
    P = map_poly_to_unit_cell(BC.B[i], BC.ll, BC.ur)
    Lx,Ly,Lz = [1.0*abs(BC.ur[a]-BC.ll[a]) for a in range(3)]
    
    # build FT terms table, using FT(-k) = FT(k)^* to skip half of calculations
    lmax = 5 # maximum l_i to sum
    FT = numpy.zeros((2*lmax+1,2*lmax+1,2*lmax+1),dtype=complex)
    p = numpy.zeros((lmax+1,lmax+1,lmax+1),dtype=complex)
    for lx in range(0,lmax+1):
        FPx = Fourier_transform_poly(P,0,0.5*lx)
        for ly in range(-lmax,lmax+1):
            FPy = Fourier_transform_poly(FPx,0,0.5*ly)
            for lz in range(-lmax,lmax+1):
                FT[lx,ly,lz] = Fourier_transform_poly(FPy,0,0.5*lz)
                FT[-lx,-ly,-lz] = FT[lx,ly,lz].conjugate()
                if lx>=0 and ly>=0 and lz>=0 and (lx,ly,lz) != (0,0,0):
                    q = sqrt( (lx/Lx)**2 + (ly/Ly)**2 + (lz/Lz)**2 )*pi
                    p[lx,ly,lz] = PG.p(q,w)

    s = 0
    eipihalf_phase = [1,1j,-1,-1j]
    for lx in range(0,lmax+1):
        for ly in range(0,lmax+1):
            for lz in range(0,lmax+1):
            
                nzero = int(lx==0)+int(ly==0)+int(lz==0)
                if nzero == 3:
                    continue    # exclude B0^2 term; also excluded by p[0,0,0]=0
                
                p_qw = p[abs(lx),abs(ly),abs(lz)]

                FTm = 0
                # only need half-sum over sx because of (sx,sy,sz)->(-sx,-sy,-sz) conjugation symmetry
                for sx in [1,]:
                    for sy in [-1,1]:
                        for sz in [-1,1]:
                            FTm += eipihalf_phase[(lx*sx+ly*sy+lz*sz)%4] * FT[lx*sx,ly*sy,lz*sz]
                FTm = (2.*FTm.real) / 8.
                
                s += p_qw * FTm**2 * 2**(3-nzero)
                                
                #print (lx,ly,lz), "\tFTm =",FTm, "\tp =",p_qw

    #print "Calculation time",(time.clock() - tstart)
    return s

def T2i(BC,PG):
    """1/T_2 dephasing time in inhomogeneous field [Hz]"""
    B0 = BC.avgB()[0]    # Average holding field, [Gauss]
    w0 = B0*PG.gmr        # Larmor angular frequency in holding field [Hz]
    
    # subtract off average holding field
    #BC.B[0] -= B0
    
    SB0B0 = SBiBi(BC,0,0,PG)
    return PG.gmr**2/4. * (2*SB0B0 + SBiBi(BC,1,w0,PG) + SBiBi(BC,2,w0,PG)).real    # 1/T_2 [Hz]





###########################
# test example calculations

def T2_Studies(outdir = os.environ["HOME"]+"/Desktop/T2_Studies/"):
    """Plots of T2 dependence on gradients and cell geometry"""
    
    PGlist = [PG_n(),PG_3He()]
    gstyles = [ [graph.style.line(),], [graph.style.line([style.linewidth.THick]),] ]
    
    ## linear gradient dBx/dx ##
    if 0:
        gT2 = graph.graphxy(width=16,height=12,
            #x=graph.axis.lin(title="Linear gradient $\\partial B_x/\\partial z$ [$\mu$G/cm]"),
            x=graph.axis.lin(title="Linear gradient $\\partial B_x/\\partial y$ [$\mu$G/cm]"),
            y=graph.axis.lin(title="Relaxation rate $1/T_2$ [mHz]",min=0),
            key = graph.key.key(pos="tl"))
        gT2.texrunner.set(lfs='foils17pt')

        gdat = []
        for dbxdz in unifrange(0,2e-6,15):
            gdat.append([dbxdz*1e6,])
            
            for PG in PGlist:
                BC = BCell()
                BC.ll,BC.ur = (-3.75,-20,-5.0),(3.75,20,5.0)    # field cell dimansions [cm]
                BC.B[0] += 0.030    # B0 = 30mG x
                BC.B[0] += monomial((0,1,0),dbxdz)            # dBx/dz gradient

                gdat[-1].append([dbxdz*1e6,1000*abs(T2i(BC,PG))])
                print PG.name #,"dbxdz = %g"%dbxdz,"1/T2 = %g"%gdat[-1][-1]
            
        for (nPG,PG) in enumerate(PGlist):
            pgdat = [p[1+nPG] for p in gdat]
            LF = LinearFitter(terms=[polyterm(2),])
            LF.fit(pgdat)
            
            gT2.plot(graph.data.points(pgdat,x=1,y=2,title="%s: $1/T_2 = %s$"%(PG.name,LF.toLatex(cfmt=".2g"))),gstyles[nPG])

        gT2.writePDFfile(outdir+"/T2_vs_dBxdy_linear.pdf")

    ## quadratic gradient dBx/dy ##
    if 1:
        gT2 = graph.graphxy(width=16,height=12,
            x=graph.axis.lin(title="Quadratic $\\langle (\\partial B_x/\\partial y)^2 \\rangle^{1/2}$ [$\mu$G/cm]"),
            y=graph.axis.lin(title="Relaxation rate $1/T_2$ [mHz]",min=0),
            key = graph.key.key(pos="tl"))
        gT2.texrunner.set(lfs='foils17pt')

        gdat = []
        for beta in unifrange(0,2e-6*sqrt(3.)/40.,15):
        
            gdat.append([beta,])
            
            for PG in PGlist:
                BC = BCell()
                BC.ll,BC.ur = (-3.75,-20,-5.0),(3.75,20,5.0)    # field cell dimansions [cm]
                BC.B[0] += 0.030    # B0 = 30mG x
                BC.B[0] += monomial((0,2,0),beta)

                gdat[-1].append([BC.GradBRMS(0,1)*1e6,1000*abs(T2i(BC,PG))])
                print PG.name,"beta = %g, <Grad^2> = %g"%(beta,(gdat[-1][-1][0]*1e-6)**2)
                print "\tuG/cm = %g, 1/T2 = %g mHz"%tuple(gdat[-1][-1])

        for (nPG,PG) in enumerate(PGlist):
            pgdat = [p[1+nPG] for p in gdat]
            LF = LinearFitter(terms=[polyterm(2),])
            LF.fit(pgdat)
            
            gT2.plot(graph.data.points(pgdat,x=1,y=2,title="%s: $1/T_2 = %s$"%(PG.name,LF.toLatex(cfmt=".2g"))),gstyles[nPG])

        gT2.writePDFfile(outdir+"/T2_vs_dBxdy_quadratic.pdf")


def ref_linear_T2_test():
    """Linear gradient reference test"""
    # From analytical calculation, should be approx.:
    #     neutron: 3.88E-5 Hz
    #     3He: 2.05E-4 Hz

    for PG in [PG_n(),PG_3He()]:
        BC = BCell()
        BC.ll,BC.ur = (-20,-3.75,-5.1),(20,3.75,5.1)    # field cell dimansions [cm]
        BC.B[0] += 0.030
        BC.B[0] += monomial((1,0,0),1e-7)
        
        print PG.name,"1/T2 = %g"%T2i(BC,PG)

    # neutron 1/T2 = -3.79892e-05
    # $^3$He 1/T2 = 0.000207002


def diagonal_multiterm_test():
    """Check answer is unchanged, for a heap of gradients"""
    for PG in [PG_n(),PG_3He()]:
        BC = BCell()
        BC.ll,BC.ur = (-20,-3.75,-5.1),(20,3.75,5.1)
        BC.B[0] += polynomial(3,{(0,0,0):0.030})
        for i in range(3):
            BC.B[0] += monomial(basisv(3,i),1e-7)
            for j in range(3):
                C = tuple([ int(a==i)+int(a==j) for a in range(3)])
                BC.B[0] += monomial(C,1e-7)
        print PG.name,"1/T2 = %g"%T2i(BC,PG)
    # for lmax = 5:
    # neutron 1/T2 = -0.00245394
    # $^3$He 1/T2 = 0.00604303


if __name__ == "__main__":
    
    diagonal_multiterm_test()
    exit(0)
    
    #ref_linear_T2_test()
    #exit(0)
    
    T2_Studies()
    exit(0)

