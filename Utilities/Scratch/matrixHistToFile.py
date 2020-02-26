# coding: utf-8
import ROOT
# Times bin size (0.5) and change units to pb
#binrescale = 0.0005
var = "pT_W"
binrescale = 10.*0.001
rtfile = ROOT.TFile("%s__NNLO_QCD.root" % var)
hist = rtfile.Get("canvas").GetListOfPrimitives().FindObject("%s__NNLO_QCD" % var)
outfile = ROOT.TFile("WplusToMuNu-MATRIX_%s.root" % var, "recreate")
outfile.mkdir("wpmunu_matrix")
outfile.cd("wpmunu_matrix")
new_hist = hist.Clone("%s_lhe_mp" % (var if var != "pT_W" else "ptW"))
new_hist.Scale(binrescale) 
new_hist.Write()
