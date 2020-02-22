# coding: utf-8
import ROOT
rtfile = ROOT.TFile("/afs/cern.ch/user/k/kelong/WplusJToENu-powheg-NNLOPS-1.root")
canvas = rtfile.Get("c")
hist = canvas.GetListOfPrimitives().FindObject("WplusJToENu-powheg-NNLOPSDYNNLO_mur1_muf1_3D.top_px")
new_canvas = ROOT.TCanvas("canvas", "canvas")
outfile = ROOT.TFile("WplusJToMuNu-powheg-NNLOPSDYNNLO_mur1_muf1.root", "recreate")
outfile.mkdir("wpmunu_dynnlo")
outfile.cd("wpmunu_dynnlo")
new_hist = hist.Clone("yW_lhe_mp")
new_hist.Scale(1.)
new_hist.Write()
