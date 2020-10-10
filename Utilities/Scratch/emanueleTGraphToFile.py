# coding: utf-8
import ROOT
import array
rtfile = ROOT.TFile("/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/wplus_rap_fewz.root")
graph = rtfile.Get("fewz")
new_canvas = ROOT.TCanvas("canvas", "canvas")
outfile = ROOT.TFile("test.root", "recreate")
outfile.mkdir("wpmunu_fewz")
outfile.cd("wpmunu_fewz")
hist = ROOT.TH1D("hist", "hist", 40, -5, 5)
for i in range(graph.GetN()):
    x = array.array('d', [0])
    y = array.array('d', [0])
    graph.GetPoint(i, x, y)
    hist.Fill(x[0], y[0])
    hist.SetBinError(hist.FindBin(x[0]), graph.GetErrorY(i))
new_hist = hist.Clone("yW_lhe_mp")
new_hist.Write()
