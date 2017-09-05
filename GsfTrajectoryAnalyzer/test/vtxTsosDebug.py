import sys
from math import sqrt,exp,log,pi
from PyMultiGaussianState1D import *
import ROOT

varAxisTitles = [ "q/p [GeV^{-1}]", "local dx/dz", "local dy/dz", "local dx [cm]", "local dy [cm]" ]

gobjects = [ ]

def drawMixture(pad,title,mgs):
    g = mgs.graph()
    g.SetLineColor(4)
    g.SetLineWidth(2)
    gobjects.append(g)

    x1,x2 = mgs.xminmax()
    mode = mgs.mode()
    y2 = mgs.mode()[1]
    hframe = pad.DrawFrame(x1-0.05*(x2-x1),y2/10000.,x2+0.05*(x2-x1),1.05*y2)
    hframe.GetXaxis().SetTitle(title)
    g.Draw("C")

    graphs = mgs.singleGraphs()
    gobjects.extend(graphs)
    for g in graphs:
        g.SetLineColor(4)
        g.SetLineStyle(2)
        g.Draw("C")
    pad.SetLogy(1)
    pad.SetGridx(1)
    pad.SetGridy(1)
    pad.Update()


if __name__=="__main__":


    cnv = ROOT.TCanvas("cnv","cnv",900,900)
    cnv.Divide(3,3)
    icnv = 0

    vtxMode = None
    coords = [ ]
    ival = None
    for il,l in enumerate(open(sys.argv[1])):

        if vtxMode==None and l.startswith("vtxTSOS mode"):
            fields = l[:-1].split()
            coord = None
            if "from" in fields[2]:
                vtxMode = fields[3][:-1]
            else:
                vtxMode = "cartesian"
            print "vtxMode = ",vtxMode
            continue

        if vtxMode!=None and ival==None:
#            assert l.startswith("printMultiState1D")
            if not l.startswith("printMultiState1D"):
                continue
            fields = l[:-1].split()
            coords.append(fields[1])
            mgs = PyMultiGaussianState1D()
            nc = int(fields[3])
            wgts = [ ]
            means = [ ]
            sigmas = [ ]
            ival = 0
            print "coord = ",coords[-1]
            continue

        if vtxMode!=None and ival!=None:
            print "Values",ival
            fields = l[:-1].strip().split()
            assert len(fields)==nc+1
            if ival==0:
                wgts = [ float(x) for x in fields[1:] ]
            elif ival==1:
                means = [ float(x) for x in fields[1:] ]
            elif ival==2:
                sigmas = [ float(x) for x in fields[1:] ]
            ival += 1 
            if ival>2:
                mgs.setStates(wgts,means,sigmas)
                icnv += 1
                cnv.cd(icnv)
                drawMixture(ROOT.gPad,coords[-1],mgs)
                ival = None
                if len(coords)==3:
                    vtxMode = None
                    coords = [ ]
            if icnv>=9:
                cnv.Update()
                raw_input("Enter")
                icnv = 0
                gobjects = [ ]



