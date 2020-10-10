#include "Analysis/SelectorTools/interface/FakeRateSelector.h"
#include <TStyle.h>

void FakeRateSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    if (!passesLeptonVeto)
        return;
    if (l1Pt < 25 || l2Pt < 15)
        return;
    if (ZMass > 101.1876 || ZMass < 81.1876)
        return;
    if (type1_pfMETEt > 25)
        return;
    if (l3MtToMET > 30)
        return;
    if (!tightZLeptons())
        return;

    if (!IsGenMatched3l())
        return;

    float pt_fillval = l3Pt;
    float eta_fillval = std::abs(l3Eta);

    float loose_weight = weight;
    if (channel_ == eee || channel_ == emm) {
        loose_weight /= eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    }
    else if (channel_ == eem || channel_ == mmm) {
        loose_weight /= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    }
    passingLoose2D_->Fill(pt_fillval, eta_fillval, loose_weight);
    passingLoose1DPt_->Fill(pt_fillval, loose_weight);
    passingLoose1DEta_->Fill(eta_fillval, loose_weight);
    if (lepton3IsTight()) {
        passingTight2D_->Fill(pt_fillval, eta_fillval, weight);
        passingTight1DPt_->Fill(pt_fillval, weight);
        passingTight1DEta_->Fill(eta_fillval, weight);
    }
}

void FakeRateSelector::SetupNewDirectory() {
    WZSelectorBase::SetupNewDirectory();

    const int nvarbins = 3;
    double variable_pt_bins[nvarbins+1] = {10, 20, 30, FR_MAX_PT_};
    AddObject<TH2D>(passingTight2D_, ("passingTight2D_"+channelName_).c_str(), "#eta; p_{T} [GeV]", nvarbins, variable_pt_bins, 3, 0, 2.5);
    AddObject<TH1D>(passingTight1DPt_, ("passingTight1DPt_"+channelName_).c_str(), "Tight leptons; p_{T} [GeV]", nvarbins, variable_pt_bins);
    AddObject<TH1D>(passingTight1DEta_, ("passingTight1DEta_"+channelName_).c_str(), "Tight leptons; #eta", 3, 0, 2.5);
    
    AddObject<TH2D>(passingLoose2D_, ("passingLoose2D_"+channelName_).c_str(), "#eta; p_{T} [GeV]", nvarbins, variable_pt_bins, 3, 0, 2.5);
    AddObject<TH1D>(passingLoose1DPt_, ("passingLoose1DPt_"+channelName_).c_str(), "Loose leptons; p_{T} [GeV]", nvarbins, variable_pt_bins);
    AddObject<TH1D>(passingLoose1DEta_, ("passingLoose1DEta_"+channelName_).c_str(), "Loose leptons; #eta", 3, 0, 2.5);
}
