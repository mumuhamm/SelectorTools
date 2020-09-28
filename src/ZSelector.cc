#include "Analysis/SelectorTools/interface/ZSelector.h"
#include <TStyle.h>
#include <regex>

void ZSelector::Init(TTree *tree)
{
    allChannels_ = {{ee, "ee"}, {mm, "mm"}};
    // Add CutFlow for Unknown to understand when channels aren't categorized
    histMap1D_[{"CutFlow", Unknown, Central}] = {};
    hists1D_ = {"CutFlow", "ZMass", "ZEta", "yZ", "ZPt", "ptl1", "etal1", "ptl2", "etal2",
		"ptj1", "ptj2", "ptj3", "etaj1", "etaj2", "etaj3", "phij1", "phij2", "phij3", "nJets",
		"MET",};

    b.SetTree(tree);
    SelectorBase::Init(tree);
    
    singleLepton_ = false;
    if (!isMC_ && name_.find("Single") != std::string::npos)
        singleLepton_ = true;
}

void ZSelector::SetScaleFactors() {
    pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
    if (pileupSF_ == nullptr ) 
        std::invalid_argument("Must pass pileup weights SF");
    eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronTightIdSF");
    if (eIdSF_ == nullptr ) 
        std::invalid_argument("Must pass electron ID SF");
    eGsfSF_ = (ScaleFactor *) GetInputList()->FindObject("electronGsfSF");
    if (eGsfSF_ == nullptr ) 
        std::invalid_argument("Must pass electron GSF SF");
    mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonTightIdSF");
    if (mIdSF_ == nullptr ) 
        std::invalid_argument("Must pass muon ID SF");
    mIsoSF_ = (ScaleFactor *) GetInputList()->FindObject("muonIsoSF");
    if (mIsoSF_ == nullptr ) 
        std::invalid_argument("Must pass muon Iso SF");

    prefireEff_ = (TEfficiency*) GetInputList()->FindObject("prefireEfficiencyMap");
    if (prefireEff_ == nullptr ) 
        std::invalid_argument("Must pass prefiring efficiency map");
}

void ZSelector::SetBranchesUWVV() {
    if (isMC_){
        fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
        fChain->SetBranchAddress("nTruePU", &nTruePU, &b_nTruePU);
    }
    fChain->SetBranchAddress("Mass", &ZMass, &b_ZMass);

    if (channel_ == ee) {
        fChain->SetBranchAddress("e1ZZTightIDNoVtx", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("e2ZZTightIDNoVtx", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("e1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("e2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("e1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("e2Eta", &l2Eta, &b_l2Eta);
    }
    else if (channel_ == mm) {
        fChain->SetBranchAddress("m1ZZTightID", &l1IsTight, &b_l1IsTight);
        fChain->SetBranchAddress("m2ZZTightID", &l2IsTight, &b_l2IsTight);
        fChain->SetBranchAddress("m1Pt", &l1Pt, &b_l1Pt);
        fChain->SetBranchAddress("m2Pt", &l2Pt, &b_l2Pt);
        fChain->SetBranchAddress("m1Eta", &l1Eta, &b_l1Eta);
        fChain->SetBranchAddress("m2Eta", &l2Eta, &b_l2Eta);
    }

    fChain->SetBranchAddress("type1_pfMETEt", &MET, &b_MET);
    fChain->SetBranchAddress("type1_pfMETPhi", &type1_pfMETPhi, &b_type1_pfMETPhi);
    fChain->SetBranchAddress("nZZTightElec", &nCBVIDTightElec, &b_nCBVIDTightElec);
    //fChain->SetBranchAddress("nCBVIDHLTSafeElec", &nCBVIDHLTSafeElec, &b_nCBVIDHLTSafeElec);
    fChain->SetBranchAddress("nZZTightMu", &nTightIdMuon, &b_nTightIdMuon);
    //fChain->SetBranchAddress("nWZMediumMuon", &nWZMediumMuon, &b_nWZMediumMuon);
}

void ZSelector::SetBranchesNanoAOD() {
    b.CleanUp();
    b.SetBranch("nElectron", nElectron);
    b.SetBranch("nMuon", nMuon);
    b.SetBranch("nJet", nJet);
    b.SetBranch("Electron_pt", Electron_pt);
    b.SetBranch("Electron_eta", Electron_eta);
    b.SetBranch("Electron_phi", Electron_phi);
    b.SetBranch("Electron_mass", Electron_mass);
    b.SetBranch("Muon_pt", Muon_pt);
    b.SetBranch("Muon_eta", Muon_eta);
    b.SetBranch("Muon_phi", Muon_phi);
    b.SetBranch("Muon_mass", Muon_mass);
    b.SetBranch("Jet_pt", Jet_pt);
    b.SetBranch("Jet_eta", Jet_eta);
    b.SetBranch("Jet_phi", Jet_phi);
    b.SetBranch("Jet_mass", Jet_mass);
    b.SetBranch("Electron_cutBased", Electron_cutBased);
    b.SetBranch("Muon_tightId", Muon_tightId);
    b.SetBranch("Muon_mediumId", Muon_mediumId);
    b.SetBranch("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
    b.SetBranch("MET_pt", MET);
    b.SetBranch("MET_phi", type1_pfMETPhi);
    b.SetBranch("Electron_charge", Electron_charge);
    b.SetBranch("Muon_charge", Muon_charge);
    //b.SetBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", Dimuon_Trigger);
    b.SetBranch("HLT_IsoMu24", SingleMuon_Trigger);
    //b.SetBranch("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", Dielectron_Trigger);
    b.SetBranch("HLT_Ele27_WPTight_Gsf", SingleElectron_Trigger);
    if (isMC_) {
        b.SetBranch("genWeight", genWeight);
        b.SetBranch("Pileup_nPU", numPU);
    }
}

void ZSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    weight = 1;
    b.SetEntry(entry);
    
    // You could resize the array per event to match the number number of muons/electrons read but
    // it doesn't seem worth it to me. Just set to a high enough value
    if (nElectron > N_KEEP_MU_E_ || nMuon > N_KEEP_MU_E_) {
        std::string message = "Found more electrons or muons than max read number.\n    Found ";
        message += std::to_string(nElectron);
        message += " electrons.\n    Found ";
        message += std::to_string(nMuon);
        message += " Muons\n  --> Max read number was ";
        message += std::to_string(N_KEEP_MU_E_);
        message += "\nExiting because this can cause problems. Increase N_KEEP_MU_E_ to avoid this error.\n";
        throw std::domain_error(message);
    }

    if (nJet > N_KEEP_JET_) {
        std::string message = "Found more jets than max read number.\n    Found ";
        message += std::to_string(nJet);
        message += " particles\n  --> Max read number was ";
        message += std::to_string(N_KEEP_JET_);
        message += "\nExiting because this can cause problems. Increase N_KEEP_JET_ to avoid this error.\n";
        throw std::domain_error(message);
    }

    ZMass = 0;
    l1Pt = 0;
    l2Pt = 0;
    l1Eta = 0;
    l2Eta = 0;
    l1Phi = 0;
    l2Phi = 0;
    l1Mass = 0;
    l2Mass = 0;

    // cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
    nCBVIDTightElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 4);
    nCBVIDVetoElec = std::count(Electron_cutBased, Electron_cutBased+nElectron, 1);
    nTightIdMuon = 0;
    nMediumIdMuon = 0;

    for (size_t i = 0; i < nMuon; i++) {
        //nMediumIdMuon += (Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15);
        nMediumIdMuon += (Muon_mediumId[i] && Muon_pfRelIso04_all[i] < 0.15);
        nTightIdMuon += (Muon_tightId[i] && Muon_pfRelIso04_all[i] < 0.15);
    }

    channel_ = channelMap_[channelName_];
    std::vector<size_t> goodIndices = {};

    if (nTightIdMuon >= 2) {
        channel_ = mm;
        channelName_ = "mm";
        if (!(Muon_mediumId[0] && Muon_pfRelIso04_all[0] < 0.15
	      && Muon_mediumId[1] && Muon_pfRelIso04_all[1] < 0.15)) {
            for (size_t i = 0; i < nMuon; i++) {
                if (Muon_mediumId[i])
                    goodIndices.push_back(i);
            }
            if (goodIndices.size() < 2) {
                return;
            }
        }
        else 
            goodIndices = {0, 1};
        if (Muon_charge[goodIndices[0]] != Muon_charge[goodIndices[1]]) {
            l1Pt = Muon_pt[goodIndices[0]];
            l2Pt = Muon_pt[goodIndices[1]];
            l1Eta = Muon_eta[goodIndices[0]];
            l2Eta = Muon_eta[goodIndices[1]];
            l1Phi = Muon_phi[goodIndices[0]];
            l2Phi = Muon_phi[goodIndices[1]];
            l1Mass = Muon_mass[goodIndices[0]];
            l2Mass = Muon_mass[goodIndices[1]];
            l1IsTight = (Muon_mediumId[goodIndices[0]] && (Muon_pfRelIso04_all[goodIndices[0]] < 0.15));
            l2IsTight = (Muon_mediumId[goodIndices[1]] && (Muon_pfRelIso04_all[goodIndices[1]] < 0.15));
        }
    }
    else if (nCBVIDTightElec >= 2) {
        channel_ = ee;
        channelName_ = "ee";
        if (!(Electron_cutBased[0] == 4 && Electron_cutBased[1] == 4)) {
            for (size_t i = 0; i < nElectron; i++) {
                if (Electron_cutBased[i] == 4)
                    goodIndices.push_back(i);
            }
            if (goodIndices.size() < 2)
                return;
        }
        else 
            goodIndices = {0, 1};
        if (Electron_charge[goodIndices[0]] != Electron_charge[goodIndices[1]]) {
            l1Pt = Electron_pt[goodIndices[0]];
            l2Pt = Electron_pt[goodIndices[1]];
            l1Eta = Electron_eta[goodIndices[0]];
            l2Eta = Electron_eta[goodIndices[1]];
            l1Phi = Electron_phi[goodIndices[0]];
            l2Phi = Electron_phi[goodIndices[1]];
            l1Mass = Electron_mass[goodIndices[0]];
            l2Mass = Electron_mass[goodIndices[1]];
            l1IsTight = (Electron_cutBased[goodIndices[0]] == 4);
            l2IsTight = (Electron_cutBased[goodIndices[1]] == 4);
        }
    }
    else {
        channel_ = Unknown;
        channelName_ = "Unknown";
    }

    leptons.clear();
    jets.clear();

    auto lep1 = LorentzVector();
    lep1.SetPt(l1Pt);
    lep1.SetEta(l1Eta);
    lep1.SetPhi(l1Phi);
    lep1.SetM(l1Mass);
    auto lep2 = LorentzVector();
    lep2.SetPt(l2Pt);
    lep2.SetEta(l2Eta);
    lep2.SetPhi(l2Phi);
    lep2.SetM(l2Mass);

    leptons.push_back(lep1);
    leptons.push_back(lep2);
    for (size_t i = 0; i < nJet; i++) {
        LorentzVector jet;
        jet.SetPt(Jet_pt[i]);
        jet.SetEta(Jet_eta[i]);
        jet.SetPhi(Jet_phi[i]);
        jet.SetM(Jet_mass[i]);
        if (jet.pt() > 30 && !helpers::overlapsCollection(jet, leptons, 0.4, leptons.size()))
            jets.push_back(jet);
    } 
    SetComposite();

    if (isMC_) {
        weight = genWeight;
        ApplyScaleFactors();
    }
    else {
        //TODO: add MET filters
    }


    //if (!singleLepton_)
    //    //passesTrigger = (Dimuon_Trigger || SingleMuon_Trigger ||
    //    //        Dielectron_Trigger || SingleElectron_Trigger);
    //    passesTrigger = (Dimuon_Trigger || SingleMuon_Trigger);
    //else
    //    passesTrigger = (!Dimuon_Trigger && SingleMuon_Trigger);
    //    //passesTrigger = ((!Dimuon_Trigger && SingleMuon_Trigger) ||
    //    //        (!Dielectron_Trigger && SingleElectron_Trigger));
    passesTrigger = SingleMuon_Trigger || SingleElectron_Trigger;

    passesLeptonVeto = (nCBVIDTightElec == 2 || nTightIdMuon == 2);
}

void ZSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation){ 
    weight = 1;
    b_l1Pt->GetEntry(entry);
    b_l2Pt->GetEntry(entry);
    b_l1Eta->GetEntry(entry);
    b_l2Eta->GetEntry(entry);

    if (isMC_) {
        b_genWeight->GetEntry(entry);
        weight = genWeight;
    }
    b_ZMass->GetEntry(entry);
    b_l1IsTight->GetEntry(entry);
    b_l2IsTight->GetEntry(entry);
    b_MET->GetEntry(entry);
    b_nCBVIDTightElec->GetEntry(entry);
    b_nTightIdMuon ->GetEntry(entry);
    
    passesTrigger = true;
    passesLeptonVeto = (nCBVIDTightElec == 2 || nTightIdMuon == 2);
}

void ZSelector::ApplyScaleFactors() {
    if (channel_ == ee) {
	if (eIdSF_ != nullptr) {
	    weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
	if (eGsfSF_ != nullptr) {
	    weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
    }
    else if (channel_ == mm) {
	if (mIdSF_ != nullptr) {
	    weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
	if (mIsoSF_ != nullptr) {
	    weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
	    weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
	}
    }
    if (pileupSF_ != nullptr) {
	weight *= pileupSF_->Evaluate1D(numPU);
    }
}

void ZSelector::SetComposite() {
    if (l1Pt == 0. || l2Pt == 0.) {
        return;
    }
    auto system = leptons.at(0)+leptons.at(1);
    ZMass = system.mass();
    ZPt = system.pt();
    ZEta = system.eta();
    Zy = system.Rapidity();
}

// Meant to be a wrapper for the tight ID just in case it changes
// To be a function of multiple variables
bool ZSelector::zlep1IsTight() {
    return l1IsTight; 
}

bool ZSelector::zlep2IsTight() {
    return l2IsTight; 
}

bool ZSelector::tightZLeptons() {
    return zlep1IsTight() && zlep2IsTight(); 
}

void ZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    int step = 0;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);


    if (channel_ != mm && channel_ != ee) 
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    if (!passesTrigger)
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    if (channel_ == ee && (std::abs(l1Eta) > 2.4 || std::abs(l2Eta) > 2.4 ))
	    return;
    else if (channel_ == mm && (std::abs(l1Eta) > 2.5 || std::abs(l2Eta) > 2.5 ))
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    if (l1Pt < 25 || l2Pt < 25)
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    if (ZMass > 106.1876 || ZMass < 76.1876)
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    //if (MET > 25)
    //    return;
    //SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    if (!tightZLeptons())
        return;
    SafeHistFill(histMap1D_, "CutFlow", channel_, variation.first, step++, weight);

    SafeHistFill(histMap1D_, "ZMass", channel_, variation.first, ZMass, weight);
    SafeHistFill(histMap1D_, "ptl1", channel_, variation.first, l1Pt, weight);
    SafeHistFill(histMap1D_, "ptl2", channel_, variation.first, l2Pt, weight);
    SafeHistFill(histMap1D_, "etal1", channel_, variation.first, l2Eta, weight);
    SafeHistFill(histMap1D_, "etal2", channel_, variation.first, l2Eta, weight);
    SafeHistFill(histMap1D_, "ZEta", channel_, variation.first, ZEta, weight);
    SafeHistFill(histMap1D_, "yZ", channel_, variation.first, Zy, weight);
    SafeHistFill(histMap1D_, "ZPt", channel_, variation.first, ZPt, weight);
    SafeHistFill(histMap1D_, "nJets", channel_, variation.first, jets.size(), weight);
    SafeHistFill(histMap1D_, "MET", channel_, variation.first, MET, weight);
    for (size_t i = 1; i <= 3; i++) {
        if (jets.size() >= i ) {
            const auto& jet = jets.at(i-1);
            SafeHistFill(histMap1D_, ("ptj"+std::to_string(i)).c_str(), channel_, variation.first, jet.pt(), weight);
            SafeHistFill(histMap1D_, ("etaj"+std::to_string(i)).c_str(), channel_, variation.first, jet.eta(), weight);
            SafeHistFill(histMap1D_, ("phij"+std::to_string(i)).c_str(), channel_, variation.first, jet.phi(), weight);
        }  
    }
}

void ZSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    InitializeHistogramsFromConfig();
}
