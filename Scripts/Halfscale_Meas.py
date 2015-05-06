#!/usr/bin/python

# Half-scale coil measurements data routines

from polynomial import *
from LinFitter import *
from StudyPlotter import *

def LF_from_poly(P):
    """Generate a linear fitter from a polynomial"""
    return LinearFitter(terms=P.getTerms(True))

def poly_from_LF(LF):
    """Re-assemble a linear fitter into a polynomial"""
    P = polynomial(LF.terms[0].N)
    for (n,t) in enumerate(LF.terms):
        P += t*LF.coeffs[n]
    return P

class ScanDataPoint:
    """One scan data point from LabView field mapper file"""
    def __init__(self,l):
    
        ###########################################################
        #    File line format:
        #
        #    0        1    2        3        4        5        6
        #    Month    Day    Year    Hour    Minute    Second    Time_s
        #
        #    7    8    9    10        11                12        13
        #    x    y    z    Current    CurrentError    1Axis    1AxisError
        #
        #    14    15        16    17        18    19
        #    Bx    BxError    By    ByError    Bz    BzError
        
        d = [float(x) for x in l.split()]
        self.abstime = d[0:6]
        self.reltime = d[6]
        self.pos = tuple([x*0.1 for x in d[7:10]])
        self.extB = d[12:14]
        self.B = (-d[14],d[16],-d[18])
        self.dB = (d[15],d[17],d[19])

    def __repr__(self):
        return "<ScanDataPt x=(%.1f, %.1f, %.1f)cm"%self.pos +"\tB=(%.3f, %.3f, %.3f)mG>"%self.B


def axis_point_groups(plist,a):
    pg = {}
    for p in plist:
        pg.setdefault(p.pos[a],[]).append(p)
    return pg

class ScanFile:

    def __init__(self,basepath,fname):
        self.basepath = basepath
        self.pts = [ ScanDataPoint(l) for l in open(basepath+"/"+fname).readlines()[1:] if len(l)>20 ]
        self.ll,self.ur = self.cell_range()
        print "Loaded",len(self.pts),"points"

    def cell_range(self):
        """Determine data points position cell range"""
        ll = [min([p.pos[a] for p in self.pts]) for a in range(3)]
        ur = [max([p.pos[a] for p in self.pts]) for a in range(3)]
        return tuple(ll), tuple(ur)
        
    def avgB(self):
        """Average field over all points"""
        return tuple([sum([p.B[a] for p in self.pts])/len(self.pts) for a in range(3)])
        
    def fit_pts(self):
        """Fit points and return a polynomial BCell"""
        
        self.LF = [LF_from_poly(lowTriangTerms(3,5)) for a in range(3)]
        
        self.BC = BCell()
        self.BC.basepath = self.basepath
        self.BC.ll,self.BC.ur = self.ll,self.ur
        
        self.Bfit = [None,None,None]
        for a in range(3):
            self.LF[a].fit([ (p.pos, p.B[a], p.dB[a]) for p in self.pts],cols=(0,1,2),errorbarWeights=True)
            print "B_%s fit RMS deviation ="%axes[a],self.LF[a].rmsDeviation()
            self.BC.B[a] = poly_from_LF(self.LF[a])

        return self.BC

    def plotCellProjection(self, xi, xj, B0=None, B0scale=None):
        
        # default to residuals plot
        if B0 is None:
            B0 = self.avgB()
        if B0scale is None:
            B0scale = sqrt(sum([x**2 for x in B0]))*0.005
            
        """Plot projection of cell fields in xi-xj plane"""
        g=graph.graphxy(width=self.ur[xi]-self.ll[xi]+1, height=self.ur[xj]-self.ll[xj]+1,
            x=graph.axis.lin(title="$%s$ position [cm]"%axes[xi], min=self.ll[xi]-0.5, max=self.ur[xi]+0.5),
            y=graph.axis.lin(title="$%s$ position [cm]"%axes[xj], min=self.ll[xj]-0.5, max=self.ur[xj]+0.5))
        g.texrunner.set(lfs='foils17pt')

        AP = ArrowPlotter(g)
        pgs = axis_point_groups(self.pts,otherAxis(xi,xj))
        pgcols = rainbowDict(pgs)
        
        for pg in pgs:
            gdat = [ (p.pos[xi],p.pos[xj], (p.B[xi]-B0[xi])/B0scale, (p.B[xj]-B0[xj])/B0scale) for p in pgs[pg]]
            asty = graph.style.arrow(lineattrs=[style.linewidth.Thick,pgcols[pg]])
            AP.plot_arrowdat(gdat,asty,offset=False)
    
        g.writePDFfile(self.basepath + "/Cell_%s-%s.pdf"%(axes[xi],axes[xj]))

    def plotFitFields(self,x0,xk, xi=None, xj=None):

        #istyles = [[rgb.red],[rgb.green],[rgb.blue]]
        #jstyles = [[style.linestyle.dotted,style.linewidth.THick],[style.linestyle.solid,style.linewidth.THick],[style.linestyle.dashed,style.linewidth.Thick]]
        
        if xi is None:
            xi = (xk+1)%3
        if xj is None:
            xj = otherAxis(xk,xi)
        
        g=graph.graphxy(width=24,height=16,
            x=graph.axis.lin(title="$%s$ position [cm]"%axes[xk]),
            y=graph.axis.lin(title="$B_{%s}$ [mG]"%axes[x0]),
            key = graph.key.key(pos="bc",columns=3))
        g.texrunner.set(lfs='foils17pt')
        
        pgsi = axis_point_groups(self.pts,xi)
        icols = rainbowDict(pgsi)
        for pgi in pgsi:
            pgsj = axis_point_groups(pgsi[pgi],xj)
            for pgj in pgsj:
                gdat = [(p.pos[xk],p.B[x0],p.dB[x0]) for p in pgsj[pgj]]
                g.plot(graph.data.points(gdat,x=1,y=2,dy=3,title=None),
                        [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[icols[pgi]]),graph.style.errorbar(errorbarattrs=[icols[pgi]])])
                        
                fpts = []
                for x in unifrange(self.ll[xk],self.ur[xk],100):
                    fpts.append([0,0,0])
                    fpts[-1][xi] = pgi
                    fpts[-1][xj] = pgj
                    fpts[-1][xk] = x
                g.plot(graph.data.points([(x[xk],self.LF[x0](x)) for x in fpts],x=1,y=2,title=None),
                        [graph.style.line(lineattrs=[icols[pgi]])])
                
        g.writePDFfile(self.basepath + "/Field_B%s_%s.pdf"%(axes[x0],axes[xk]))

    def plot_noise_v_time(self):
    
        g=graph.graphxy(width=24,height=16,
            x=graph.axis.lin(title="scan time [m]",min=0,max=max([p.reltime/60. for p in self.pts])),
            y=graph.axis.lin(title="$\\Delta B$ [$\\mu$G]"),
            key = graph.key.key(pos="bc",columns=3))
        g.texrunner.set(lfs='foils17pt')
        
        acols = [rgb.red,rgb.green,rgb.blue]
        B0 = self.avgB()
        for a in range(3):
            if 0:
                gdat = [ (p.reltime,1000*p.dB[a]) for p in self.pts]
                g.plot(graph.data.points(gdat,x=1,y=2,title="$B_%s$"%axes[a]),
                            [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[acols[a]])])
            #gdat = [ (p.reltime,p.B[a]-B0[a],p.dB[a]) for p in self.pts]
            gdat = [ (p.reltime/60.,1000*(p.B[a]-self.LF[a](p.pos)),1000*p.dB[a]) for p in self.pts]
            g.plot(graph.data.points(gdat,x=1,y=2,dy=3,title="$B_%s$"%axes[a]),
                        [graph.style.symbol(symbol.circle,size=0.2,symbolattrs=[acols[a]]),graph.style.errorbar(errorbarattrs=[acols[a]])])
            
        
        g.writePDFfile(self.basepath + "/B_spread_v_time.pdf")


if __name__=="__main__":

    F = ScanFile("/Users/michael/Desktop/EDM_Field_Maps/20140306_partial_cooldown","EDMFieldMapping2014-03-06_15.31.txt")
    BC = F.fit_pts()

    # display average gradients
    print "B0 =",F.avgB()
    print "gradient = (%.2f, %.2f, %.2f) uG/cm"%tuple([1000*x for x in BC.GradBAvg()])

    if 0:
        for x0 in range(3):
            for xk in range(3):
                F.plotFitFields(x0,xk)
            
    F.plot_noise_v_time()
                
    if 0:
        F.plotCellProjection(0,2) # x-z
        F.plotCellProjection(0,1) # x-y
        F.plotCellProjection(1,2) # y-z



# x,y speed limit 30, 90 for z; units in mm
# meas time, ~10s
# 1s wait time after every move; dwell time



