#!/usr/bin/env python
import ROOT
import argparse
import os
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import time
from python import UserInput,HistTools,ConfigureJobs
import array

def getHist(rtfile, dataset, histname, xsec, rebin):
    hist = rtfile.Get("/".join([dataset, histname]))
    sumw = rtfile.Get("/".join([dataset, "sumweights"]))
    if sumw:
        hist.Scale(xsec/sumw.Integral())

    if rebin:
        hist = HistTools.rebinHist(hist, histname+"rebin", rebin)
    print "Dataset", dataset, "Integral is", hist.Integral()
    return hist

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, help="Name of output scale factor")
parser.add_argument("--hist_names", type=str, nargs=2, help="Name of histograms (numerator denominator)")
parser.add_argument("--append", action='store_true', help="Add to file (by default overwrite file)")
parser.add_argument("-a", "--analysis", type=str, help="Analysis name")
parser.add_argument("-n", "--numerator", type=str, help="Name of numerator dataset")
parser.add_argument("-d", "--denominator", type=str, help="Name of denominator dataset")
parser.add_argument("inputFile", type=str, help="file with histograms")
parser.add_argument("-r", "--rebin", 
                    type=str, default=None, help="Rebin array: "
                    "values (bin edges) separated by commas.")
args = parser.parse_args()

args.rebin = array.array('d', UserInput.getRebin(args.rebin))

output_file = os.environ["CMSSW_BASE"]+'/src/Analysis/SelectorTools/data/%s_scaleFactors.root' % args.analysis
fScales = ROOT.TFile(output_file, 'recreate' if not args.append else 'update')

rtfile = ROOT.TFile.Open(args.inputFile)

xsecs  = ConfigureJobs.getListOfFilesWithXSec([args.numerator, args.denominator], 
                                                os.path.expanduser("~/work/"))

histnum = getHist(rtfile, args.numerator, args.hist_names[0], xsecs[args.numerator], args.rebin)
histdenom = getHist(rtfile, args.denominator, args.hist_names[1], xsecs[args.denominator], args.rebin)

hist = histnum.Clone(args.name)
hist.Divide(histdenom)

sf = ROOT.ScaleFactor(args.name, args.name)
sf.Set1DHist(hist, 0, 0)
fScales.cd()
sf.Write()
