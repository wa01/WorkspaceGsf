import sys

import ROOT

tf = ROOT.TFile(sys.argv[1])
analyzerDir = tf.Get("trajectoryAnalyzer")

for ev in analyzerDir.Get("GsfTree"):
    print ""
    print ev.run,ev.lumi,ev.evt,ev.itraj,ev.itmf,ev.itmr
    print ev.fwPredLPar[0],ev.fwPredLErr[0]
    ev.fwPredGPos.Print()
