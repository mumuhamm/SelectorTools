#include "Analysis/VVAnalysis/interface/WGenSelector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <TStyle.h>
#include <regex>
#include <cmath>

void WGenSelector::Init(TTree *tree)
{
    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    allChannels_ = {{ep, "ep"}, {en, "en"}, {mp, "mp"}, {mn, "mn"}};
    hists1D_ = {"CutFlow", "mWmet", "yWmet", "ptWmet", "mW", "yW", "ptW", "mTtrue", "mTmet",
        "ptl", "etal", "phil", "ptnu", "etanu", "phinu", "MET", "MET_phi",
        "ptj1", "ptj2", "etaj1", "etaj2", "nJets",
        "dRlgamma_maxptassoc", "dRlgamma_minassoc", "ptg_closeassoc", "ptg_maxassoc", "nGammaAssoc", 
        "ptgmax_assoc", "ptgmax_assoc",
    };
    doSystematics_ = true;
    systematics_ = {
        {BareLeptons, "barelep"},
        {BornParticles, "born"},
        {LHEParticles, "lhe"},
    };
    systHists_ = hists1D_;

    weighthists1D_ = {"CutFlow", "yW", "ptW", "mTtrue", "ptl", "etal", "phil", "ptnu", "etanu", };

    nLeptons_ = 1;
    doNeutrinos_ = true;
    doPhotons_ = true;
    doFiducial_ = true;
    doTheoryVars_ = true;
    NanoGenSelectorBase::Init(tree);
}

void WGenSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);

    if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 11) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = en;
            channelName_ = "en";
        }
        else {
            channel_ = ep;
            channelName_ = "ep";
        }
    }
    else if (leptons.size() > 0 && std::abs(leptons.at(0).pdgId()) == 13) {
        if (leptons.at(0).pdgId() > 0) {
            channel_ = mn;
            channelName_ = "mn";
        }
        else {
            channel_ = mp;
            channelName_ = "mp";
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
        return;
    }
}

void WGenSelector::SetComposite() {
    if (leptons.size() == 0) {
        wCandMet = LorentzVector();
        wCand = LorentzVector();
        return;
    }
    else if (neutrinos.size() == 0) {
        wCand = LorentzVector();
        return;
    }
    auto lepP4 = leptons.at(0).polarP4();
    auto compareByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() < b.pt(); };
    auto mt = [] (LorentzVector& l, LorentzVector& v) {
        return std::sqrt(2*l.pt()*v.pt()*(1 - cos(l.phi() - v.phi())));
    };

    auto nup = std::max_element(neutrinos.begin(), neutrinos.end(), compareByPt);
    if (neutrinos.size()) {
        nu = neutrinos.size() > 0 ? nup->polarP4() : LorentzVector();
        wCandMet = lepP4 + genMet;
        wCand = neutrinos.size() > 0 ? lepP4 + nu : LorentzVector();
        mTtrue = mt(lepP4, nu);
    }
    mTmet = mt(lepP4, genMet);
}

void WGenSelector::FillHistograms(Long64_t entry, SystPair variation) { 
    std::string lepType = "";
    FillHistogramsByName(entry, lepType, variation);
}

void WGenSelector::FillHistogramsByName(Long64_t entry, std::string& toAppend, SystPair variation) { 
    int step = 0;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (channel_ != mn && channel_ != en && channel_ != mp && channel_ != ep) 
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (leptons.size() < nLeptons_)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (neutrinos.size() < nLeptons_)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    auto& lep = leptons.at(0);
    if (doFiducial_ && std::abs(lep.eta()) > 2.5)
        return;
    SafeHistFill(histMap1D_, concatenateNames("CutFlow", toAppend), channel_, variation.first, step++, weight);

    if (doFiducial_ && lep.pt() < 25)
        return;

    SafeHistFill(histMap1D_, concatenateNames("mW", toAppend), channel_, variation.first, wCand.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yW", toAppend), channel_, variation.first, wCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptW", toAppend), channel_, variation.first, wCand.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("mTtrue", toAppend), channel_, variation.first, mTtrue, weight);
    SafeHistFill(histMap1D_, concatenateNames("mTmet", toAppend), channel_, variation.first, mTmet, weight);
    SafeHistFill(histMap1D_, concatenateNames("mWmet", toAppend), channel_, variation.first, wCandMet.mass(), weight);
    SafeHistFill(histMap1D_, concatenateNames("yWmet", toAppend), channel_, variation.first, wCandMet.Rapidity(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptWmet", toAppend), channel_, variation.first, wCandMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET", toAppend), channel_, variation.first, genMet.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("MET_phi", toAppend), channel_, variation.first, genMet.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptl", toAppend), channel_, variation.first, lep.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("etal", toAppend), channel_, variation.first, lep.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phil", toAppend), channel_, variation.first, lep.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("ptnu", toAppend), channel_, variation.first, nu.pt(), weight);
    SafeHistFill(histMap1D_, concatenateNames("etanu", toAppend), channel_, variation.first, nu.eta(), weight);
    SafeHistFill(histMap1D_, concatenateNames("phinu", toAppend), channel_, variation.first, nu.phi(), weight);
    SafeHistFill(histMap1D_, concatenateNames("nJets", toAppend), channel_, variation.first, jets.size(), weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, concatenateNames(("ptj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("etaj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("phij"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.phi(), weight);
        }  
    }

    if (variation.first == Central && doTheoryVars_) {
        size_t maxEntry = *nLHEScaleWeight+*nLHEPdfWeight;
        if (doMC2H_ == true)
            maxEntry += N_MC2HESSIAN_WEIGHTS_;

        for (size_t i = 0; i < maxEntry; i++) {
            //float thweight = (i < *nLHEScaleWeight) ? LHEScaleWeight[i] : ( i < (*nLHEPdfWeight+*nLHEScaleWeight) ? LHEPdfWeight[i-*nLHEScaleWeight] 
            //        : LHEHessianPdfWeight[i-*nLHEScaleWeight-*nLHEPdfWeight]);
            float thweight = i < *nLHEScaleWeight ? LHEScaleWeight[i] : LHEPdfWeight[i-*nLHEScaleWeight];
            if (i > *nLHEScaleWeight && i < *nLHEScaleWeight+*nLHEPdfWeight)
                mcPdfWeights_->Fill(thweight);
            else if (i > *nLHEScaleWeight+*nLHEPdfWeight)
                hesPdfWeights_->Fill(thweight);
            else 
                scaleWeights_->Fill(thweight);
            thweight *= weight;
            SafeHistFill(weighthistMap1D_, concatenateNames("mW", toAppend), channel_, variation.first, wCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yW", toAppend), channel_, variation.first, wCand.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptW", toAppend), channel_, variation.first, wCand.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("mWmet", toAppend), channel_, variation.first, wCandMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("yWmet", toAppend), channel_, variation.first, wCandMet.Rapidity(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptWmet", toAppend), channel_, variation.first, wCandMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("MET", toAppend), channel_, variation.first, genMet.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("MET_phi", toAppend), channel_, variation.first, genMet.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("ptl", toAppend), channel_, variation.first, lep.pt(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("etal", toAppend), channel_, variation.first, lep.eta(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("phil", toAppend), channel_, variation.first, lep.phi(), i, thweight);
            SafeHistFill(weighthistMap1D_, concatenateNames("nJets", toAppend), channel_, variation.first, jets.size(), i, thweight);
        }
    }

    if (variation.first == BareLeptons) {
        SafeHistFill(histMap1D_, "nGammaAssoc", channel_, variation.first, photons.size(), weight);

        auto compareByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() < b.pt(); };
        auto compareByDRLead = [lep] (const reco::GenParticle& a, const reco::GenParticle& b) {
            return reco::deltaR(a, lep) < reco::deltaR(b, lep);
        };

        auto gclose = std::min_element(photons.begin(), photons.end(), compareByDRLead);
        auto maxPtg = std::max_element(photons.begin(), photons.end(), compareByPt);

        SafeHistFill(histMap1D_, "dRlgamma_minassoc", channel_, variation.first, photons.size() > 0 ? reco::deltaR(*gclose, lep) : 0., weight);
        SafeHistFill(histMap1D_, "dRlgamma_maxptassoc", channel_, variation.first, photons.size() > 0 ? reco::deltaR(*maxPtg, lep) : 0., weight);
        SafeHistFill(histMap1D_, "ptg_closeassoc", channel_, variation.first, photons.size() > 0 ? gclose->pt() : 0., weight);
        SafeHistFill(histMap1D_, "ptgmax_assoc", channel_, variation.first, photons.size() > 0 ? maxPtg->pt() : 0., weight);
    }
}
