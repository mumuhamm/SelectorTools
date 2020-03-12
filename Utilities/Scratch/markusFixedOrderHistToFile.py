# coding: utf-8
import ROOT
#name = "WplusJToENu"
name = "ZJToEE"
chan = "mm"
outname = "DYm50"
filename = "/eos/user/k/kelong/DYNNLOFiles/%s-powheg-NNLOPS-1.root" % name
print filename
rtfile = ROOT.TFile(filename)
canvas = rtfile.Get("c")
hist = canvas.GetListOfPrimitives().FindObject(name+"-powheg-NNLOPSDYNNLO_mur1_muf1_3D.top_px")
new_canvas = ROOT.TCanvas("canvas", "canvas")
outfile = ROOT.TFile("/eos/user/k/kelong/HistFiles/ZGen/NoSelection/%s-powheg-NNLOPSDYNNLO_mur1_muf1.root" % name, "recreate")
outdir = "%s_dynnlo" % outname
outfile.mkdir(outdir)
outfile.cd(outdir)
new_hist = hist.Clone("%s_lhe_%s" % ("yW" if "WJ" in name else "yZ", chan))
new_hist.Scale(1.)
new_hist.Write()
