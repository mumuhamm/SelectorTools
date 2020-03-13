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
    hists2D_ = {"etal_ptl_2D"};
    doSystematics_ = true;
    systematics_ = {
        {mWShift100MeVUp, "mWShift100MeVUp"},
        {mWShift50MeVUp, "mWShift50MeVUp"},
        {mWShift25MeVUp, "mWShift25MeVUp"},
        {mWShift20MeVUp, "mWShift20MeVUp"},
        {mWShift10MeVUp, "mWShift10MeVUp"},
        {mWShift100MeVDown, "mWShift100MeVDown"},
        {mWShift50MeVDown, "mWShift50MeVDown"},
        {mWShift25MeVDown, "mWShift25MeVDown"},
        {mWShift20MeVDown, "mWShift20MeVDown"},
        {mWShift10MeVUp, "mWShift10MeVUp"},
        {BareLeptons, "barelep"},
        {BornParticles, "born"},
        {LHEParticles, "lhe"},
        {muonScaleUp, "CMS_scale_mUp"},
        {muonScaleDown, "CMS_scale_mDown"},
    };
    systHists_ = hists1D_;
    systHists2D_ = hists2D_;

    weighthists1D_ = {"CutFlow", "yW", "ptW", "mTtrue", "ptl", "etal", "phil", "ptnu", "etanu", };
    weighthists2D_ = hists2D_;
    //weightsysts1D_ = {LHEParticles};
    //weightsysts2D_ = {LHEParticles};

    refWeight = 1;
    nLeptons_ = 1;
    doNeutrinos_ = true;
    doPhotons_ = true;
    
    // Chose by MC sample
    if (name_.find("nnlops") != std::string::npos) {
        MW_GEN_ = 80398.0;
        GAMMAW_GEN_ = 2088.720;
    }
    else if (name_.find("minnlo") != std::string::npos) {
        MW_GEN_ = 80398.0;
        GAMMAW_GEN_ = 2141.0;
    }
    else {
        MW_GEN_ = 80419.;
        GAMMAW_GEN_ = 2050;
    }

    NanoGenSelectorBase::Init(tree);
}

void WGenSelector::LoadBranchesNanoAOD(Long64_t entry, SystPair variation) { 
    NanoGenSelectorBase::LoadBranchesNanoAOD(entry, variation);
    if (variation.first == Central)
        cenWeight = weight;
    else if (variation.first == muonScaleUp && leptons.size() >= nLeptons_) {
        auto& l = leptons.at(0);
        leptons[0] = makeGenParticle(l.pdgId(), l.status(), l.pt()*1.0001, l.eta(), l.phi(), l.mass());
        SetComposite();
    }
    else if (variation.first == muonScaleDown && leptons.size() >= nLeptons_) {
        auto& l = leptons.at(0);
        leptons[0] = makeGenParticle(l.pdgId(), l.status(), l.pt()*0.9999, l.eta(), l.phi(), l.mass());
        SetComposite();
    }
    else if (variation.first == LHEParticles)
        mWlhe = wCand.mass()*1000.;
    else if (variation.first == mWShift10MeVUp)
        weight = cenWeight*breitWignerWeight(10.);
    else if (variation.first == mWShift10MeVDown)
        weight = cenWeight*breitWignerWeight(-10.);
    else if (variation.first == mWShift20MeVUp)
        weight = cenWeight*breitWignerWeight(20.);
    else if (variation.first == mWShift20MeVDown)
        weight = cenWeight*breitWignerWeight(-20.);
    else if (variation.first == mWShift25MeVUp)
        weight = cenWeight*breitWignerWeight(25.);
    else if (variation.first == mWShift25MeVDown)
        weight = cenWeight*breitWignerWeight(-25.);
    else if (variation.first == mWShift50MeVUp)
        weight = cenWeight*breitWignerWeight(50.);
    else if (variation.first == mWShift50MeVDown)
        weight = cenWeight*breitWignerWeight(-50.);
    else if (variation.first == mWShift100MeVUp)
        weight = cenWeight*breitWignerWeight(100.);
    else if (variation.first == mWShift100MeVDown)
        weight = cenWeight*breitWignerWeight(-100.);

    if (variation.first == LHEParticles) {
        ptVlhe = wCand.pt();
    }

    if (variation.first == ptV0to3 && (ptVlhe < 0. || ptVlhe > 3.))
        leptons.clear();
    else if (variation.first == ptV3to5 && (ptVlhe < 3. || ptVlhe > 5.))
        leptons.clear();
    else if (variation.first == ptV5to7 && (ptVlhe < 5. || ptVlhe > 7.))
        leptons.clear();
    else if (variation.first == ptV7to9 && (ptVlhe < 7. || ptVlhe > 9.))
        leptons.clear();
    else if (variation.first == ptV9to12 && (ptVlhe < 9. || ptVlhe > 12.))
        leptons.clear();
    else if (variation.first == ptV12to15 && (ptVlhe < 12. || ptVlhe > 15.))
        leptons.clear();
    else if (variation.first == ptV15to20 && (ptVlhe < 15. || ptVlhe > 20.))
        leptons.clear();
    else if (variation.first == ptV20to27 && (ptVlhe < 20. || ptVlhe > 27.))
        leptons.clear();
    else if (variation.first == ptV27to40 && (ptVlhe < 27. || ptVlhe > 40.))
        leptons.clear();
    else if (variation.first == ptV40toInf && ptVlhe > 40. )
        leptons.clear();

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

double WGenSelector::breitWignerWeight(double offset) {

    double targetMass = MW_GEN_ + offset;
    double s_hat = mWlhe*mWlhe;
    double offshell = s_hat - MW_GEN_*MW_GEN_;
    double offshellOffset = s_hat - targetMass*targetMass;
    double weight = (offshell*offshell + GAMMAW_GEN_*GAMMAW_GEN_*MW_GEN_*MW_GEN_)/
            (offshellOffset*offshellOffset + GAMMAW_GEN_*GAMMAW_GEN_*targetMass*targetMass);
    return weight;
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
    mcWeights_->Fill(weight/refWeight);
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
    SafeHistFill(histMap2D_, concatenateNames("etal_ptl_2D", toAppend), channel_, variation.first, lep.eta(), lep.pt(), weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, concatenateNames(("ptj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("etaj"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, concatenateNames(("phij"+std::to_string(i)).c_str(), toAppend), channel_, variation.first, jet.phi(), weight);
        }  
    }

    if (std::find(theoryVarSysts_.begin(), theoryVarSysts_.end(), variation.first) != theoryVarSysts_.end()) {
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
            if (centralWeightIndex_ != -1)
                thweight /= LHEScaleWeight.At(centralWeightIndex_);
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
            SafeHistFill(weighthistMap2D_, concatenateNames("etal_ptl_2D", toAppend), channel_, variation.first, lep.eta(), lep.pt(), i, thweight);
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
