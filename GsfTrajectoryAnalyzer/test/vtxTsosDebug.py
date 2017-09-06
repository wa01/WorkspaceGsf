import sys
from math import sqrt,exp,log,pi
from PyMultiGaussianState1D import *
import ROOT

varAxisTitles = [ "q/p [GeV^{-1}]", "local dx/dz", "local dy/dz", "local dx [cm]", "local dy [cm]" ]

gobjects = [ ]

def drawMixture(pad,title,mgs,xmin=None,xmax=None,superimpose=False):
    g = mgs.graph()
    color = 4 if not superimpose else 2
    g.SetLineColor(color)
    g.SetLineWidth(2)
    gobjects.append(g)

    x1,x2 = mgs.xminmax()
    mode = mgs.mode()
    y2 = mgs.mode()[1]
    if xmin==None:
        xmin = x1-0.05*(x2-x1)
    if xmax==None:
        xmax = x2+0.05*(x2-x1)
    if not superimpose:
        hframe = pad.DrawFrame(xmin,y2/5000.,xmax,2.*y2)
        hframe.GetXaxis().SetTitle(title)
    g.Draw("C")

    graphs = mgs.singleGraphs()
    gobjects.extend(graphs)
    for g in graphs:
        g.SetLineColor(color)
        g.SetLineStyle(2)
        g.Draw("C")
    pad.SetLogy(1)
    pad.SetGridx(1)
    pad.SetGridy(1)
    pad.Update()

    line = "Mean : {0:10.2e}#pm{1:9.2e}  Mode : {2:10.2e}".format(mgs.mean(),mgs.sigma(),mode[0])
    if superimpose:
        latex = ROOT.TLatex(0.1,0.91,line)
        latex.SetTextColor(2)
    else:
        latex = ROOT.TLatex(0.1,0.95,line)
        latex.SetTextColor(4)
    latex.SetNDC(1)
#        pave = ROOT.TPaveText(0.1,0.80,0.5,0.9,"NDC NB")
    latex.SetTextSize(0.04)
    latex.SetTextAlign(11)
    latex.Draw()
    gobjects.append(latex)

class DrawSingleTrack:

    def __init__(self,fname,index=None,otherCanvas=None,event=(None,None,None)):
        self.fname = fname
        self.index = index
        self.otherCanvas = otherCanvas
        self.canvasTitle = None

        self.event = event
        self.runLumiEvt = None
        self.itrack = 0
        self.canvas = None
        self.icnv = 0

        self.vtxMode = None
        self.coords = [ ]
        self.ival = None

        self.gobjects = [ ]
        

    def drawTrack(self):

        if self.otherCanvas==None:
            self.canvas = ROOT.TCanvas("cnv","cnv",900,900)
            self.canvas.Divide(3,3)
        else:
            self.canvas = self.otherCanvas
        self.icnv = 0

        self.vtxMode = None
        self.coords = [ ]
        self.ival = None
        showEvent = True
        for il,l in enumerate(open(self.fname)):

            if self.icnv==0:
                px, py, pz = None, None, None

            if self.otherCanvas==None and l.startswith("Begin processing"):
                fields = l[:-1].split()
                if fields[5]=="Run" and fields[7]=="Event" and fields[9]=="LumiSection":
                    self.canvasTitle = "Run "+fields[6][:-1]+" Lumi "+fields[10]+" Evt "+fields[8][:-1]
                    self.runLumiEvt = ( int(fields[6][:-1]), int(fields[10]), int(fields[8][:-1]) )
                    showEvent = ( self.event[0]==None or self.event[0]==self.runLumiEvt[0] ) and \
                                ( self.event[1]==None or self.event[1]==self.runLumiEvt[1] ) and \
                                ( self.event[2]==None or self.event[2]==self.runLumiEvt[2] )
                continue

            if self.vtxMode==None and l.startswith("vtxTSOS mode"):
                fields = l[:-1].split()
                coord = None
                if "from" in fields[2]:
                    self.vtxMode = fields[3][:-1]
                else:
                    self.vtxMode = "cartesian"
                print "vtxMode = ",self.vtxMode
                continue

            if self.vtxMode!=None and self.ival==None:
    #            assert l.startswith("printMultiState1D")
                if not l.startswith("printMultiState1D"):
                    continue
                fields = l[:-1].split()
                self.coords.append(fields[1])
                mgs = PyMultiGaussianState1D()
                nc = int(fields[3])
                wgts = [ ]
                means = [ ]
                sigmas = [ ]
                self.ival = 0
                # print "coord = ",self.coords[-1]
                continue

            if self.vtxMode!=None and self.ival!=None:
                # print "Values",self.ival
                fields = l[:-1].strip().split()
                assert len(fields)==nc+1
                if self.ival==0:
                    wgts = [ float(x) for x in fields[1:] ]
                elif self.ival==1:
                    means = [ float(x) for x in fields[1:] ]
                elif self.ival==2:
                    sigmas = [ float(x) for x in fields[1:] ]
                self.ival += 1 
                if self.ival>2:
                    mgs.setStates(wgts,means,sigmas)
                    self.icnv += 1
                    self.canvas.cd(self.icnv)
                    if ( self.index==None or self.itrack==self.index ) and showEvent:
                        if self.coords[-1] in [ 'px', 'py', 'pz', 'p' ]:
                            drawMixture(ROOT.gPad,self.coords[-1],mgs,xmin=-200.,xmax=200., \
                                        superimpose=(self.otherCanvas!=None))
                            if self.coords[-1]=='px':
                                px = mgs.mode()[0]
#                                print "px = ",px
                            elif self.coords[-1]=='py':
                                py = mgs.mode()[0]
#                                print "py = ",py
                            elif self.coords[-1]=='pz':
                                pz = mgs.mode()[0]
#                                print "pz = ",pz
                        elif self.coords[-1]=='qp':
                            drawMixture(ROOT.gPad,self.coords[-1],mgs,xmin=-1./2.5,xmax=1./2.5, \
                                        superimpose=(self.otherCanvas!=None))
                        else:
                            drawMixture(ROOT.gPad,self.coords[-1],mgs, \
                                        superimpose=(self.otherCanvas!=None))
                        ROOT.gPad.Update()
                    self.ival = None
                    if len(self.coords)==3:
                        self.vtxMode = None
                        self.coords = [ ]
                if self.icnv>=9:
                    if self.canvasTitle!=None:
                        title = "Track "+str(self.itrack)
                        if px!=None and py!=None and pz!=None:
                            p3 = ROOT.TVector3(px,py,pz)
                            title += " ; Pt {0:4.1f} Eta {1:5.2f}".format(p3.Perp(),p3.Eta())
                        title += " ; "+self.canvasTitle
                        self.canvas.SetTitle(title)
                    if ( self.index==None or self.itrack==self.index ) and showEvent:
                        self.canvas.Update()
                        raw_input("Enter")
                    if self.itrack==self.index:
                        return self.canvas
                    else:
                        self.itrack += 1
                    self.icnv = 0
                    self.gobjects = [ ]



if __name__=="__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--itrack', '-i', type=int, default=None)
    parser.add_argument('--run', '-r', type=int, default=None)
    parser.add_argument('--lumi', '-l', type=int, default=None)
    parser.add_argument('--evt', '-e', type=int, default=None)
    parser.add_argument('files',nargs='+',default=None)
    args = parser.parse_args()
    assert len(args.files)<=2 and (len(args.files)==1 or args.itrack!=None)

    drawClass = DrawSingleTrack(fname=args.files[0],index=args.itrack,event=(args.run,args.lumi,args.evt))
    cnv = drawClass.drawTrack()

    if len(args.files)>1:    
        drawClass2 = DrawSingleTrack(fname=args.files[1],index=args.itrack,otherCanvas=cnv)
        cnv2 = drawClass2.drawTrack()

