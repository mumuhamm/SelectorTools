#include "Analysis/VVAnalysis/interface/ZGenSelector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TStyle.h>
#include <regex>

void ZGenSelector::Init(TTree *tree)
{
    allChannels_ = {{ee, "ee"}, {mm, "mm"}};
    // Add CutFlow for Unknown to understand when channels aren't categorized
    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    std::vector<std::string> basehists1D = {"CutFlow", "ZMass", "yZ", "ptZ", "ptl1", "etal1", "phil1", "ptl2", "etal2", "phil2", 
        "ptj1", "ptj2", "ptj3", "etaj1", "etaj2", "etaj3", "phij1", "phij2", "phij3", "nJets",
        "MET", "HT",};
    hists1D_ = basehists1D;
    std::vector<std::string> partonicChans = {"uu_dd", "uubar_ddbar", "ug_dg", "ubarg_dbarg", "gg", "other"};
    for (auto& chan : partonicChans) {
        for (auto& hist : basehists1D)
            hists1D_.push_back(chan + "_" + hist);
    }

    doSystematics_ = true;
    systematics_ = {
        {BareLeptons, "barelep"},
        {BornParticles, "born"},
        {LHEParticles, "lhe"},
    };
    systHists_ = hists1D_;

    weighthists1D_ = {"CutFlow", "ZMass", "yZ", "ptZ", "ptl1", "etal1", "phil1", "ptl2", "etal2", "phil2", 
        "ptj1", "ptj2", "ptj3", "etaj1", "etaj2", "etaj3", "phij1", "phij2", "phij3", "nJets",
        "MET", "HT", };
    nLeptons_ = 2;

    NanoGenSelectorBase::Init(tree);
    
}

void ZGenSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);

    if (leptons.size() < 2) {
        channel_ = Unknown;
        channelName_ = "Unknown";
        return;
    }
    if (leptons.at(0).pdgId() + leptons.at(1).pdgId() == 0) {
        if (std::abs(leptons.at(0).pdgId()) == 11) {
            channel_ = ee;
            channelName_ = "ee";
        }
        else if (std::abs(leptons.at(0).pdgId()) == 13) {
            channel_ = mm;
            channelName_ = "mm";
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
    }
}

void ZGenSelector::SetComposite() {
    if (leptons.size() >= 2)
        zCand = leptons.at(0).polarP4() + leptons.at(1).polarP4();
}

void ZGenSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    int step = 0;
    int failStep = 0;

    step++;
    if (channel_ != mm && channel_ != ee) 
        failStep = step;

    auto lep1 = leptons.size() > 1 ? leptons.at(0) : reco::GenParticle();
    auto lep2 = leptons.size() > 1 ? leptons.at(1) : reco::GenParticle();

    step++;
    if (zCand.mass() < 50.)
        failStep = step;
    step++;
    if (lep1.pt() < 25. || lep2.pt() < 25.)
        failStep = step;
    step++;
    if (std::abs(lep1.eta()) > 2.5 || std::abs(lep2.eta()) > 2.5)
        failStep = step;
    step++;
    if (zCand.mass() < 60. || zCand.mass() > 120.)
        failStep = step;

    for (int j = 0; j < (failStep == 0 ? step : failStep); j++) {
        SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, j, weight);
        for (size_t i = 0; i < *nLHEScaleWeight+*nLHEPdfWeight; i++) {
            float thweight = i < *nLHEScaleWeight ? LHEScaleWeight[i] : LHEPdfWeight[i-*nLHEScaleWeight];
            thweight *= weight;
            SafeHistFill(weighthistMap1D_, "CutFlow", channel_, variation.first, j, i, thweight);
        }
    }
    if (doFiducial_ && failStep != 0)
        return;

    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);
    SafeHistFill(histMap1D_, "ZMass", channel_, variation.first, zCand.mass(), weight);
    SafeHistFill(histMap1D_, "yZ", channel_, variation.first, zCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, "ptZ", channel_, variation.first, zCand.pt(), weight);
    SafeHistFill(histMap1D_, "ptl1", channel_, variation.first, lep1.pt(), weight);
    SafeHistFill(histMap1D_, "etal1", channel_, variation.first, lep1.eta(), weight);
    SafeHistFill(histMap1D_, "phil1", channel_, variation.first, lep1.phi(), weight);
    SafeHistFill(histMap1D_, "ptl2", channel_, variation.first, lep2.pt(), weight);
    SafeHistFill(histMap1D_, "etal2", channel_, variation.first, lep2.eta(), weight);
    SafeHistFill(histMap1D_, "phil2", channel_, variation.first, lep2.phi(), weight);
    SafeHistFill(histMap1D_, "nJets", channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap1D_, "MET", channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, "HT", channel_, variation.first, ht, weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, ("ptj"+std::to_string(i)).c_str(), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, ("etaj"+std::to_string(i)).c_str(), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, ("phij"+std::to_string(i)).c_str(), channel_, variation.first, jet.phi(), weight);
        }  
    }
    if (variation.first == Central) {
        for (size_t i = 0; i < *nLHEScaleWeight+*nLHEPdfWeight; i++) {
            float thweight = i < *nLHEScaleWeight ? LHEScaleWeight[i] : LHEPdfWeight[i-*nLHEScaleWeight];
            thweight *= weight;
            SafeHistFill(weighthistMap1D_, "ZMass", channel_, variation.first, zCand.mass(), i, thweight);
            SafeHistFill(weighthistMap1D_, "yZ", channel_, variation.first, zCand.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptZ", channel_, variation.first, zCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptl1", channel_, variation.first, lep1.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "etal1", channel_, variation.first, lep1.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, "phil1", channel_, variation.first, lep1.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, "ptl2", channel_, variation.first, lep2.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, "etal2", channel_, variation.first, lep2.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, "phil2", channel_, variation.first, lep2.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, "nJets", channel_, variation.first, jets.size(), i, thweight);
            SafeHistFill(weighthistMap1D_, "MET", channel_, variation.first, genMet.pt(), i, thweight);
        }
    }
    // Should check how slow this is

    std::string partonicChan = "other";
    if ((*Generator_id1 == 1 && *Generator_id2 == 1) || (*Generator_id1 == 2 && *Generator_id2 == 2))
        partonicChan = "uu_dd";
    else if ((*Generator_id1 == 1 && *Generator_id2 == -1) || (*Generator_id1 == 2 && *Generator_id2 == -2))
        partonicChan = "uubar_ddbar";
    else if (*Generator_id1 == 21 && *Generator_id2 == 21)
        partonicChan = "gg";
    else if ((*Generator_id1 == 1 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == 1) || 
                (*Generator_id1 == 2 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == 2))
        partonicChan = "ug_dg";
    else if ((*Generator_id1 == -1 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == -1) || 
                (*Generator_id1 == -2 && *Generator_id2 == 21) || (*Generator_id1 == 21 && *Generator_id2 == -2))
        partonicChan = "ubarg_dbarg";

    SafeHistFill(histMap1D_, (partonicChan+"_ZMass").c_str(), channel_, variation.first, zCand.mass(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_yZ").c_str(), channel_, variation.first, zCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptZ").c_str(), channel_, variation.first, zCand.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptl1").c_str(), channel_, variation.first, lep1.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_etal1").c_str(), channel_, variation.first, lep1.eta(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_phil1").c_str(), channel_, variation.first, lep1.phi(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_ptl2").c_str(), channel_, variation.first, lep2.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_etal2").c_str(), channel_, variation.first, lep2.eta(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_phil2").c_str(), channel_, variation.first, lep2.phi(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_nJets").c_str(), channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_MET").c_str(), channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, (partonicChan+"_HT").c_str(), channel_, variation.first, ht, weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, (partonicChan+"_ptj"+std::to_string(i)).c_str(), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, (partonicChan+"_etaj"+std::to_string(i)).c_str(), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, (partonicChan+"_phij"+std::to_string(i)).c_str(), channel_, variation.first, jet.phi(), weight);
        }  
    }
}

