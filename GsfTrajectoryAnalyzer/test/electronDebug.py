import sys
from TrajectoryStates import *
import ROOT

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
    assert ev.ic==-1
    tm = Measurement(ev.itmf)
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
        print ncFwPred,ncBwPred,ncUpd,ev.ic,ev.wgtFwPred,ev.wgtBwPred,ev.wgtUpd
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
    tm.setBwPred(bwPred)
    tm.setUpd(upd)
    return tm


def readTrajectory(iev):
    ev = iev.next(keep=True)
    assert ev.itmf==0
    traj = Trajectory(ev.run,ev.lumi,ev.evt,ev.itraj)

#    endFlg = False
    for i in range(-ev.itmr+1):
        tm = readMeasurement(iev)
#        if endFlg:
#            assert i==-ev.itmr
#            break

    return traj

if __name__=="__main__":


        
    tf = ROOT.TFile(sys.argv[1])
    analyzerDir = tf.Get("trajectoryAnalyzer")
    gsfTree = analyzerDir.Get("GsfTree")

    iev = MyIterator(gsfTree)
#    endFlg = False
    while True:
        readTrajectory(iev)
        break


