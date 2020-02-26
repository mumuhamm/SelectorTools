# coding: utf-8
import ROOT
rtfile = ROOT.TFile("yW__NNLO_QCD.root")
hist = rtfile.Get("canvas").GetListOfPrimitives().FindObject("yW__NNLO_QCD")
outfile = ROOT.TFile("WplusToMuNu-MATRIX.root", "recreate")
outfile.mkdir("wpmunu_matrix")
outfile.cd("wpmunu_matrix")
new_hist = hist.Clone("yW_lhe_mp")
# Times bin size (0.5) and change units to pb
print new_hist.Integral()
new_hist.Scale(0.0005) 
print new_hist.Integral()
new_hist.Write()
