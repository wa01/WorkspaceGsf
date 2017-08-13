import sys
from TrajectoryStates import *
import ROOT

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
        assert self.Upd==None
        self.Upd = state

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

def readMeasurement(ev):
    assert ev.ic==-1
    tm = Measurement(ev.imf)
    fwPred = MultiTState(parameters=ev.fwPredLPar,errors=ev.fwPredLErr)
    bwPred = MultiTState(parameters=ev.bwPredLPar,errors=ev.bwPredLErr)
    upd = MultiTState(parameters=ev.updLPar,errors=updLErr)
    ncFwPred = ev.ncFwPred
    ncBwPred = ev.ncBwPred
    ncUpd = ev.ncUpd
    ncMax = max(ncFwPred,ncBwPred,ncUpd)
    endFlg = False
    for i in range(ncMax):
        if i<ncFwPred:
            assert ev.wgtFwPred>=0
            fwPred.addComponent(SingleTState(ev.fwPredLPar,ev.fwPredLErr,ev.wgtFwPred))
        if i<ncBwPred:
            assert ev.wgtBwPred>=0
            bwPred.addComponent(SingleTState(ev.bwPredLPar,ev.bwPredLErr,ev.wgtBwPred))
        if i<ncUpd:
            assert ev.wgtUpd>=0
            upd.addComponent(SingleTState(ev.updLPar,ev.updLErr,ev.wgtUpd))
        try:
            next(ev)
        except StopIteration:
            assert i==(ncMax-1)
            endFlg = True
            break
    tm.setFwPred(fwPred)
    tm.setBwPred(bwPred)
    tm.setUpd(upd)
    return ( endFlg, tm )


def readTrajectory(ev):
    assert ev.imf==0
    traj = Trajectory(ev.run,ev.lumi,ev.evt,ev.itraj)

    for i in range(-ev.imr+1):
        endFlg, tm = readMeasurement(ev)
        if endFlg:
            assert i==-ev.imr



        
