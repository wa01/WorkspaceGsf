# import ROOT in batch mode
import sys
from math import sqrt,exp,log,pi

varAxisTitles = [ "q/p [GeV^{-1}]", "local dx/dz", "local dy/dz", "local dx [cm]", "local dy [cm]" ]

#oldargv = sys.argv[:]
#sys.argv = [ '-b-' ]
import ROOT
#ROOT.gROOT.SetBatch(True)
#sys.argv = oldargv

OneOverSqrt2Pi = 1./sqrt(2*pi)

def evalSingleGauss(x,weight,mean,sigma):
    dx = (x-mean)/sigma
    if abs(dx)>20.:
        return 0.
    return weight*OneOverSqrt2Pi/sigma*exp(-0.5*dx*dx)

def evalMultiGauss(x,weights,means,sigmas):
    result = 0.
    for i in range(len(weights)):
        result += evalSingleGauss(x,weights[i],means[i],sigmas[i])
    return result

class GaussianMixture1D:

    def __init__(self,weights,means,sigmas):
        self.nc_ = len(weights)
        self.weights_ = weights
        self.means_ = means
        assert len(means)==self.nc_ and len(sigmas)==self.nc_
        self.sigmas_ = sigmas
        self.graph_ = None
        self.singleGraphs_ = None
        self.mode_ = None
        self.mean_ = None
        self.sigma_ = None
    
    def eval(self,x):
        return evalMultiGauss(x,self.weights_,self.means_,self.sigmas_)

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

        if self.nc_==1:
            self.mean_ = self.means_[0]
            self.sigma_ = self.sigmas_[0]
            return

        meanMean = 0.
        weightSum = 0.
        measVar1 = 0.
        measVar2 = 0.
        for i1 in range(self.nc_):
            weight = self.weights_[i1]
            weightSum += weight

            mean1 = self.means_[i1]
            meanMean += weight * mean1
            measVar1 += weight * self.sigmas_[i1]**2

            for i2 in range(i1+1,self.nc_):
                posDiff = mean1 - self.means_[i2]
                measVar2 += weight * posDiff**2

        self.mean_ = meanMean/weightSum
        self.sigma_ = sqrt(measVar1/weightSum+measVar2/weightSum/weightSum)

    def mean(self):
        if self.mean_==None:
            self.computeMean()
        return self.mean_

    def sigma(self):
        if self.sigma_==None:
            self.computeMean()
        return self.sigma_

             
    def xminmax(self,sigmaMax=5.):
        result = [ None, None ]
        for i in range(len(self.weights_)):
            xl = self.means_[i] - sigmaMax*self.sigmas_[i]
            if result[0]==None or xl<result[0]:
                result[0] = xl
            xh = self.means_[i] + sigmaMax*self.sigmas_[i]
            if result[1]==None or xh>result[1]:
                result[1] = xh
        return tuple(result)
        
    def graph(self,deltaSigma=0.25,sigmaMax=5.,recalculate=True):
        if self.graph_!=None and not recalculate:
            return self.graph_

        minsig = min(self.sigmas_)
        points = [ ]
        for ic in range(self.nc_):
            sx = -sigmaMax
            while sx<(1.+deltaSigma/2.)*sigmaMax:
                x = self.means_[ic] + sx*self.sigmas_[ic]
                points.append( ( x, evalMultiGauss(x,self.weights_,self.means_,self.sigmas_) ) )
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
        for ic in range(self.nc_):
            g = ROOT.TGraph()
            np = 0
            sx = -sigmaMax
            while sx<(1.+deltaSigma/2.)*sigmaMax:
                x = self.means_[ic] + sx*self.sigmas_[ic]
                g.SetPoint(np,x,evalSingleGauss(x,self.weights_[ic],self.means_[ic],self.sigmas_[ic]))
                np += 1
                sx += deltaSigma
            self.singleGraphs_.append(g)

        return self.singleGraphs_
        
# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

gsfElectronHandle, gsfElectronLabel = Handle("std::vector<reco::GsfElectron>"), ("gedGsfElectrons")

events = Events(sys.argv[1])


for iev,event in enumerate(events):

    print ""
    print ""
    print event.eventAuxiliary().run(),event.eventAuxiliary().luminosityBlock(),event.eventAuxiliary().event()
    event.getByLabel(gsfElectronLabel,gsfElectronHandle)

    if (iev%10000)==0:
        print "At event ",iev

    print ""
    latex = ROOT.TLatex()
    latex.SetTextSize(0.03)
    latex.SetNDC(1)
#    marker = ROOT.TMarker()
#    marker.SetMarkerStyle(20)
    for iEle,gsfEle in enumerate(gsfElectronHandle.product()):
        gsfTrack = gsfEle.gsfTrack()
#        if abs(gsfTrack.eta()-gsfTrack.etaMode())>0.001 or abs(gsfTrack.phi()-gsfTrack.phiMode())>0.001:
        print "electron pt,eta,phi      = {0:8.1f} {1:7.3f} {2:7.2f}".format(gsfEle.pt(),gsfEle.eta(),gsfEle.phi())
        print "gsfTrack pt,eta,phi      = {0:8.1f} {1:7.3f} {2:7.2f}".format(gsfTrack.pt(),gsfTrack.eta(),gsfTrack.phi())
        print "gsfTrack mode pt,eta,phi = {0:8.1f} {1:7.3f} {2:7.2f}".format(gsfTrack.ptMode(),gsfTrack.etaMode(),gsfTrack.phiMode())
        print "gsfTrack lostHits / validHits / nchi2 = ",gsfTrack.numberOfLostHits(),gsfTrack.numberOfValidHits(),gsfTrack.normalizedChi2()
        print " std /  corrected cluster energy = ",gsfEle.ecalEnergy(),gsfEle.correctedEcalEnergy()
        print " fbrem = ",gsfEle.fbrem()

        if len(sys.argv)>=3 and event.eventAuxiliary().event()==int(sys.argv[2]) and \
                ( len(sys.argv)==3 or iEle==int(sys.argv[3]) ):
            gsfExtra = gsfTrack.gsfExtra()
            weights = gsfExtra.innerStateWeights()
            ROOT.gROOT.cd()
            cnv = ROOT.TCanvas("cnv","cnv",600,900)
            cnv.Divide(2,3)
            graphs = [ ]
            for iv in range(5):
                if iv==0:
                    pad = cnv.cd(1)
                else:
                    pad = cnv.cd(iv+2)
                means = [ x[iv] for x in gsfExtra.innerStateLocalParameters() ]
                sigmas = [ sqrt(x[iv][iv]) for x in gsfExtra.innerStateCovariances() ]
                mixture = GaussianMixture1D(weights,means,sigmas)
                g = mixture.graph()
                graphs.append(g)
                x1,x2 = mixture.xminmax()
                mode = mixture.mode()
                hframe = pad.DrawFrame(x1-0.05*(x2-x1),mode[1]/1000000.,x2+0.05*(x2-x1),1.05*mode[1])
                hframe.GetXaxis().SetTitle(varAxisTitles[iv])
                g.Draw("C")
                for gs in mixture.singleGraphs():
                    gs.SetLineStyle(2)
                    gs.Draw("C")
                    graphs.append(gs)
#                marker.DrawMarker(mode[0],mode[1])
                latex.DrawLatex(0.65,0.85,"Mode  = {0:10.5f}".format(mode[0]))
                latex.DrawLatex(0.65,0.80,"Mean  = {0:10.5f}".format(mixture.mean()))
                latex.DrawLatex(0.65,0.75,"Sigma = {0:10.5f}".format(mixture.sigma()))
                pad.Update()
            cnv.Update()
            raw_input("Enter")
