import ROOT
from math import sqrt,exp,log,pi

class PyMultiGaussianState1D:

    def __init__(self):
        self.OneOverSqrt2Pi = 1./sqrt(2*pi)

        if not ROOT.gROOT.GetClass("MultiGaussianState1D"):
            print "Not found"
            ROOT.gROOT.ProcessLine(".L MultiGaussianState1D.cc+")
        self.singleStates = ROOT.MultiGaussianState1D.SingleState1dContainer()
        self.nc = None
        self.multiState_ = None

        self.graph_ = None
        self.singleGraphs_ = None
        self.mode_ = None
        self.mean_ = None
        self.sigma_ = None

    def addState(self,weight,mean,sigma):
        assert self.multiState_==None
        self.singleStates.push_back(ROOT.SingleGaussianState1D(mean,sigma**2,weight))

    def setStates(self,weights,means,sigmas):
        assert self.multiState_==None
        assert len(weights)==len(means) and len(weights)==len(sigmas)
        for i in range(len(weights)):
            self.addState(weights[i],means[i],sigmas[i])

    def nComponents(self):
        if self.nc==None:
            self.nc = self.singleStates.size()
            self.multiState_ = ROOT.MultiGaussianState1D(self.singleStates)
        return self.nc

    def multiState(self):
        if self.multiState_==None:
            self.nc = self.singleStates.size()
            self.multiState_ = ROOT.MultiGaussianState1D(self.singleStates)
        return self.multiState_

    def evalComponent(self,x,i):
        assert i>=0 and i<self.nComponents()
        mean = self.multiState().components()[i].mean()
        sigma = self.componentSigma(i)
        dx = (x-mean)/sigma
        if abs(dx)>20.:
            return 0.
        weight = self.multiState().components()[i].weight()
        return weight*self.OneOverSqrt2Pi/sigma*exp(-0.5*dx*dx)

    def eval(self,x):
        result = 0.
        for i in range(self.nComponents()):
            result += self.evalComponent(x,i)
        return result

    def mode(self):
        if self.mode_!=None:
            return self.mode_
        g = self.graph()
        np = g.GetN()
        xs = g.GetX()
        ys = g.GetY()
        iymax = None
        ymax = 0.
        for ip in range(np):
            if ys[ip]>ymax:
                iymax = ip
                ymax = ys[ip]
        assert iymax!=None and iymax>0 and iymax<(np-1)
        x1,y1 = xs[iymax-1],ys[iymax-1]
        x2,y2 = xs[iymax],ys[iymax]
        x3,y3 = xs[iymax+1],ys[iymax+1]
        a = ((y3-y2)/(x3-x2)-(y2-y1)/(x2-x1)) / (x3-x1)
        b = (y2-y1)/(x2-x1) - a*(x2+x1)
        c = y1 - a*x1*x1 - b*x1
        assert a<0.
        self.mode_ = ( -b/(2*a), self.eval(-b/(2*a)) )
        return self.mode_

    def computeMean(self):
        self.mean_ = self.multiState().mean()
        self.sigma_ = sqrt(self.multiState().variance())

    def mean(self):
        if self.mean_==None:
            self.computeMean()
        return self.mean_

    def sigma(self):
        if self.sigma_==None:
            self.computeMean()
        return self.sigma_
             
    def componentWeight(self,i):
        assert i>=0 and i<self.nComponents()
        return self.multiState().components()[i].weight()

    def componentMean(self,i):
        assert i>=0 and i<self.nComponents()
        return self.multiState().components()[i].mean()

    def componentSigma(self,i):
        assert i>=0 and i<self.nComponents()
        return self.multiState().components()[i].standardDeviation()

    def xminmax(self,sigmaMax=5.,min=None,max=None):
        result = [ None, None ]
        if min!=None and max!=None:
            result = [ min, max ]
        for i in range(self.nComponents()):
            xl = self.componentMean(i) - sigmaMax*self.componentSigma(i)
            if result[0]==None or xl<result[0]:
                result[0] = xl
            xh = self.componentMean(i) + sigmaMax*self.componentSigma(i)
            if result[1]==None or xh>result[1]:
                result[1] = xh
        return tuple(result)
        
    def graph(self,deltaSigma=0.25,sigmaMax=5.,recalculate=True):
        if self.graph_!=None and not recalculate:
            return self.graph_

        minsig = min([ self.componentSigma(i) for i in range(self.nComponents()) ])
        points = [ ]
        for ic in range(self.nComponents()):
            sx = -sigmaMax
            while sx<(1.+deltaSigma/2.)*sigmaMax:
                x = self.componentMean(ic) + sx*self.componentSigma(ic)
                points.append( ( x, self.eval(x) ) )
                sx += deltaSigma
        points.sort()
        self.graph_ = ROOT.TGraph()
        np = 0
        xlast = None
        for x,v in points:
            if xlast==None or (x-xlast)>0.8*deltaSigma*minsig:
                self.graph_.SetPoint(np,x,v)
                np += 1
                xlast = x
        return self.graph_

    def singleGraphs(self,deltaSigma=0.25,sigmaMax=5.,recalculate=True):
        if self.singleGraphs_!=None and not recalculate:
            return self.singleGraphs_

        self.singleGraphs_ = [ ]
        for ic in range(self.nComponents()):
            g = ROOT.TGraph()
            np = 0
            sx = -sigmaMax
            while sx<(1.+deltaSigma/2.)*sigmaMax:
                x = self.componentMean(ic) + sx*self.componentSigma(ic)
                g.SetPoint(np,x,self.evalComponent(x,ic))
                np += 1
                sx += deltaSigma
            self.singleGraphs_.append(g)

        return self.singleGraphs_
        
