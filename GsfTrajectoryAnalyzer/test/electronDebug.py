import sys
from math import sqrt,exp,log,pi
from TrajectoryStates import *
import ROOT

varAxisTitles = [ "q/p [GeV^{-1}]", "local dx/dz", "local dy/dz", "local dx [cm]", "local dy [cm]" ]

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
        if type(weights)==type([]):
            self.nc_ = len(weights)
            self.weights_ = weights
            self.means_ = means
            self.sigmas_ = sigmas
            assert len(means)==self.nc_ and len(sigmas)==self.nc_
        else:
            assert type(means)==type(weights) and type(sigmas)==type(weights)
            self.nc = 1
            self.weights_ = [ weights ]
            self.means_ = [ means ]
            self.sigmas_ = [ sigmas ]
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

             
    def xminmax(self,sigmaMax=5.,min=None,max=None):
        result = [ None, None ]
        if min!=None and max!=None:
            result = [ min, max ]
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
        

class MyIterator:
    def __init__(self,itObj):
        self.iter = iter(itObj)
        self.current = None
        self.fakeNext = False

    def next(self,keep=False):
        if not self.fakeNext:
            self.current = next(self.iter)
        self.fakeNext = False
        if keep:
            self.fakeNext = True
        return self.current

    def current(self):
        return self.current

class Measurement:
    def __init__(self,im):
        self.im = im
        self.hit = None
        self.bwPred = None
        self.fwPred = None
        self.upd = None

    def setHit(self,hit):
        assert self.hit==None
        self.hit = hit

    def setBwPred(self,state):
        assert self.bwPred==None
        self.bwPred = state

    def setFwPred(self,state):
        assert self.fwPred==None
        self.fwPred = state

    def setUpd(self,state):
        assert self.upd==None
        self.upd = state

class Trajectory:
    def __init__(self,run,lumi,evt,itraj):
        self.run = run
        self.lumi = lumi
        self.evt = evt
        self.itraj = itraj

        self.nm = 0
        self.measurements = [ ]

    def addMeasurement(self,measurement):
        self.nm += 1
        self.measurements.append(measurement)

def readMeasurement(iev):
    ev = iev.next()
#    print "Next",ev.run,ev.lumi,ev.itmf,ev.itmr,ev.ic,ev.wgtFwPred,ev.wgtBwPred,ev.wgtUpd
#    print "    ",ev.fwPredLPar[0],ev.fwPredLErr[0],ev.bwPredLPar[0],ev.bwPredLErr[0], \
#        ev.updLPar[0],ev.updLErr[0]
    assert ev.ic==-1
    tm = Measurement(ev.itmf)
    tm.setHit(SingleTState(parameters=ev.hitPar,errors=ev.hitErr))
    fwPred = MultiTState(parameters=ev.fwPredLPar,errors=ev.fwPredLErr)
    bwPred = MultiTState(parameters=ev.bwPredLPar,errors=ev.bwPredLErr)
    upd = MultiTState(parameters=ev.updLPar,errors=ev.updLErr)
    ncFwPred = ev.ncFwPred
    ncBwPred = ev.ncBwPred
    ncUpd = ev.ncUpd
    ncMax = max(ncFwPred,ncBwPred,ncUpd)
#    endFlg = False
    for i in range(ncMax):
        ev = iev.next()
#        print "Next",ev.run,ev.lumi,ev.itmf,ev.itmr,ev.ic,ev.wgtFwPred,ev.wgtBwPred,ev.wgtUpd
#        print "    ",ev.fwPredLPar[0],ev.fwPredLErr[0],ev.bwPredLPar[0],ev.bwPredLErr[0], \
#            ev.updLPar[0],ev.updLErr[0]
#        print ncFwPred,ncBwPred,ncUpd,ev.ic,ev.wgtFwPred,ev.wgtBwPred,ev.wgtUpd
        if i<ncFwPred:
            assert ev.wgtFwPred>=0
            fwPred.addComponent(ev.fwPredLPar,ev.fwPredLErr,ev.wgtFwPred)
        if i<ncBwPred:
            assert ev.wgtBwPred>=0
            bwPred.addComponent(ev.bwPredLPar,ev.bwPredLErr,ev.wgtBwPred)
        if i<ncUpd:
            assert ev.wgtUpd>=0
            upd.addComponent(ev.updLPar,ev.updLErr,ev.wgtUpd)
#        ev = None
#        try:
#            next(ev)
#        except StopIteration:
#            assert i==(ncMax-1)
#            endFlg = True
#            break
    tm.setFwPred(fwPred)
#    print tm.fwPred.components_[0].pars_[0],tm.fwPred.components_[0].errors_[0]
    tm.setBwPred(bwPred)
    tm.setUpd(upd)
    return tm


def readTrajectory(iev):
    ev = iev.next(keep=True)
#    print "Next",ev.run,ev.lumi,ev.itmf,ev.itmr,ev.ic,ev.wgtFwPred,ev.wgtBwPred,ev.wgtUpd
#    print "    ",ev.fwPredLPar[0],ev.fwPredLErr[0],ev.bwPredLPar[0],ev.bwPredLErr[0], \
#        ev.updLPar[0],ev.updLErr[0]
    assert ev.itmf==0
    traj = Trajectory(ev.run,ev.lumi,ev.evt,ev.itraj)

#    endFlg = False
    for i in range(-ev.itmr+1):
#        print "itm",i
        tm = readMeasurement(iev)
#        print tm.fwPred.components_[0].pars_[0],tm.fwPred.components_[0].errors_[0]
        traj.addMeasurement(tm)
#        print traj.measurements[0].fwPred.components_[0].pars_[0],traj.measurements[0].fwPred.components_[0].errors_[0]
#        print traj.measurements[0].bwPred.components_[0].pars_[0],traj.measurements[0].bwPred.components_[0].errors_[0]
#        if endFlg:
#            assert i==-ev.itmr
#            break

#    print traj.measurements[0].fwPred.components_[0].pars_[0],traj.measurements[0].fwPred.components_[0].errors_[0]
    return traj

if __name__=="__main__":


    iv = 0

    tf = ROOT.TFile(sys.argv[1])
    analyzerDir = tf.Get("trajectoryAnalyzer")
    gsfTree = analyzerDir.Get("GsfTree")

    iev = MyIterator(gsfTree)
#    endFlg = False
    while True:
        traj = readTrajectory(iev)
#        print traj

#        ROOT.gROOT.cd()
        cnv = ROOT.TCanvas("cnv","cnv",800,1000)
        cnv.Divide(3,int((traj.nm+0.5)/3))
        graphs = [ ]

        for im in range(traj.nm):
            print im
            tm = traj.measurements[im]
            if iv>2:
                peGauss = tm.hit[iv-3]
                hitMix = GaussianMixture1D(1.,peGauss[0],peGauss[1])
            peGauss = tm.fwPred[iv]
            fwPredMix = GaussianMixture1D(tm.fwPred.weights(),peGauss[0],peGauss[1])
            print "fwPred"
            print tm.fwPred.weights()
            print peGauss[0]
            print peGauss[1]
            print tm.fwPred.components_[0].pars_[iv],tm.fwPred.components_[0].errors_[iv]
            peGauss = tm.bwPred[iv]
            bwPredMix = GaussianMixture1D(tm.bwPred.weights(),peGauss[0],peGauss[1])
            print "bwPred"
            print tm.bwPred.weights()
            print peGauss[0]
            print peGauss[1]
            if im==0 or im==(traj.nm-1):
                peGauss = tm.upd[iv]
                updMix = GaussianMixture1D(tm.upd.weights(),peGauss[0],peGauss[1])
                print "upd"
                print tm.upd.weights()
                print peGauss[0]
                print peGauss[1]
            else:
                updMix = None
            gFw = fwPredMix.graph()
            gFw.SetLineColor(4)
            graphs.append(gFw)
            gBw = bwPredMix.graph()
            gBw.SetLineColor(2)
            graphs.append(gBw)
            if updMix!=None:
                gUpd = updMix.graph()
                gUpd.SetLineColor(1)
                graphs.append(gUpd)

            x1f,x2f = fwPredMix.xminmax()
            fwPredMode = fwPredMix.mode()
            x1b,x2b = bwPredMix.xminmax()
            if (x2f-x1f)<(x2b-x1b):
                x1,x2 = x1f,x2f
            else:
                x1,x2 = x1b,x2b
            bwPredMode = bwPredMix.mode()
            y2 = max(fwPredMode[1],bwPredMode[1])
            if updMix!=None:
#                x1,x2 = updMix.xminmax(min=x1,max=x2)
                updMode = updMix.mode()
                y2 = max(y2,updMode[1])

            cnv.cd(im+1)
            print x1-0.05*(x2-x1),y2/1000000.,x2+0.05*(x2-x1),1.05*y2
            hframe = ROOT.gPad.DrawFrame(x1-0.05*(x2-x1),y2/10000.,x2+0.05*(x2-x1),1.05*y2)
            hframe.GetXaxis().SetTitle(varAxisTitles[iv])
            gFw.Draw("C")
            gBw.Draw("C")
            if updMix!=None:
                gUpd.Draw("C")
            ROOT.gPad.SetLogy(1)
            ROOT.gPad.Update()

        cnv.Update()

        raw_input("Enter")
        break


