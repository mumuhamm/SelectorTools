import ROOT
from python import HistTools

names = [
    91.0876,
    91.0926,
    91.0976,
    91.1026,
    91.1076,
    91.1126,
    91.1176,
    91.1226,
    91.1276,
    91.1326,
    91.1376,
    91.1426,
    91.1476,
    91.1526,
    91.1576,
    91.1626,
    91.1676,
    91.1726,
    91.1776,
    91.1826,
    91.1876,
    91.1926,
    91.1976,
    91.2026,
    91.2076,
    91.2126,
    91.2176,
    91.2226,
    91.2276,
    91.2326,
    91.2376,
    91.2426,
    91.2476,
    91.2526,
    91.2576,
    91.2626,
    91.2676,
    91.2726,
    91.2776,
    91.2826,
    91.2876,
]

def getMass(hist):
    fit = hist.GetFunction("gaus")
    return fit.GetParameter(1)

rtfile = ROOT.TFile("/eos/user/k/kelong/HistFiles/ZGen/NoSelection/ZMiNNLOUpdate_EWWeights.root")
hist = rtfile.Get("DYm50_minnlo_ewweights__photos/ZMass_mm")
hist2D = rtfile.Get("DYm50_minnlo_ewweights__photos/ZMass_lheWeights_mm")

hist.Fit("gaus", "Q", "", 81.2, 101.2)
cenPeak = getMass(hist)
print "-"*80
print "Peak for central 91.1876 (peak at 91.1535) is %0.4f" % cenPeak 
print "-"*80

hists,_ = HistTools.getLHEWeightHists(hist2D, range(1,100), "", "mZshift")

for i,h in enumerate(hists[18:18+len(names)]):
    print "-"*80
    print "Results for weight %i for mass shift %0.4f" % (i, names[i])
    print "    --> offset from central should be %0.4f" % (names[i] - 91.1876)
    print "-"*80
    h.Fit("gaus", "Q", "", 80, 100)
    mass = getMass(h)
    print "Peak is at %0.4f offset from central is %0.4f" % (mass, mass - cenPeak)
