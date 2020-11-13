# Setup ScaleFactor objects to be "registered" with ROOT,
# allowing them to be called from TTree.Draw(), for example.
# Currently used for lepton scale factors and pileup weights.
#
# Modified from N. Smith, U. Wisconsin
# 

#!/usr/bin/env python
import ROOT
import argparse
import os
#ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import time

from python import sfhelpers

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tightfr_file", type=str,
        default='data/fakeRate18Apr2017-3LooseLeptons-TightMuons.root')
parser.add_argument("-m", "--medfr_file", type=str,
        default='data/fakeRate18Apr2017-3LooseLeptons-MediumMuons.root')
args = parser.parse_args()
output_file = 'data/scaleFactors.root'
fScales = ROOT.TFile(output_file, 'recreate')

# For nTruePU reweighting
pileupSF = ROOT.ScaleFactor("pileupSF", "Run2016B-H 36.8/fb Pileup profile over RunIISpring16 MC Scale Factor, x=NTruePU")
pileupFile = ROOT.TFile.Open('PileupWeights/PU_Central.root')
pileupFileUp = ROOT.TFile.Open('PileupWeights/PU_minBiasUP.root')
pileupFileDown = ROOT.TFile.Open('PileupWeights/PU_minBiasDOWN.root')
pileupSF.Set1DHist(pileupFile.Get('pileup'), pileupFileUp.Get('pileup'), pileupFileDown.Get('pileup'))
fScales.cd()
pileupSF.Write()


MuonScaleList = [["data/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root","SF"],
                 ["data/TnP_NUM_MultiIsoMedium_DENOM_MediumID_VAR_map_pt_eta.root","SF"],
                 ["data/TnP_NUM_PtError_DENOM_LooseID_VAR_map_pt_eta.root","SF"],
                 ["data/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root","SF"],
                 ["data/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root","SF"]]

ElectronScaleList = [["data/Electron_scaleFactors.root","GsfElectronToMVATightIDEmuTightIP2DSIP3D4"],
                     ["data/Electron_scaleFactors.root","MVATightElectronToConvVetoIHit0"],
                     ["data/Electron_scaleFactors.root","MVATightConvIHit0ElectronToCharge"],
                     ["data/Electron_scaleFactors.root","MVATightElectronToMultiIsoTISOEmu"],]


muonIdSF = ROOT.ScaleFactor("muonMediumIdSF", "Moriond '17 Muon Tight WP ID SF, x=abs(Eta), y=Pt, z=run number") 

sfhelper.setScaleFactorObj(MuonScaleList, muonIdSF, fScales)

ElecIdSF = ROOT.ScaleFactor("electronTightIdSF", "Moriond '17 Electron Tight WP ID SF, x=abs(Eta), y=Pt, z=run number")

sfhelper.setScaleFactorObj(ElectronScaleList, ElecIdSF, fScales)
