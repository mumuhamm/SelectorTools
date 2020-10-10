#include "Analysis/SelectorTools/interface/WZSelectorBase.h"
#include "TLorentzVector.h"
#include <TStyle.h>
#include <regex>
#include "TParameter.h"

// This is very WZ specific and should really be improved or likely removed
std::string WZSelectorBase::GetNameFromFile() {
    std::regex expr = std::regex("201[0-9]-[0-9][0-9]-[0-9][0-9]-(.*)-WZxsec2016");
    std::smatch matches;
    std::string fileName = fChain->GetTree()->GetDirectory()->GetFile()->GetName(); 

    std::regex_search(fileName, matches, expr);
    return std::string(matches.str(1));
}

void WZSelectorBase::SetScaleFactors() {
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

void WZSelectorBase::Init(TTree *tree)
{
    isVBS_ = (selection_ == VBSselection_Loose || 
        selection_ == VBSselection_Loose_Full || 
        selection_ == VBSselection_NoZeppenfeld || 
        selection_ == VBSselection_NoZeppenfeld_Full || 
        selection_ == VBSselection_Tight || 
        selection_ == VBSselection_Tight_Full || 
        selection_ == VBSBackgroundControl || 
        selection_ == VBSBackgroundControl_Full || 
        selection_ == VBSBackgroundControlATLAS || 
        selection_ == VBSBackgroundControlLoose ||
        selection_ == VBSBackgroundControlLoose_Full
        );

    allChannels_ = {{eee, "eee"}, {eem, "eem"}, {emm, "emm"}, {mmm, "mmm"}};

    if (isMC_){
        isNonpromptMC_ = false;
        isZgamma_ = false;
        if (std::find(nonprompt3l_.begin(), nonprompt3l_.end(), name_) != nonprompt3l_.end()) {
            isNonpromptMC_ = true;
        }
        else if (name_ == "zg") {
            isZgamma_ = true;
        }
    }

    b.SetTree(tree);
    SelectorBase::Init(tree);
}

void WZSelectorBase::SetBranchesUWVV() {
    b.CleanUp();
    if (isMC_){
        b.SetBranch("genWeight", genWeight);
        b.SetBranch("nTruePU", nTruePU);
    }
    else {
        b.SetBranch("Flag_duplicateMuonsPass", Flag_duplicateMuonsPass);
        b.SetBranch("Flag_badMuonsPass", Flag_badMuonsPass);
    }

    if (channel_ == eee) {
        b.SetBranch("e1IsCBVIDTight", l1IsTight);
        b.SetBranch("e2IsCBVIDTight", l2IsTight);
        b.SetBranch("e3IsCBVIDTight", l3IsTight);
        b.SetBranch("e1PVDXY", l1PVDXY);
        b.SetBranch("e2PVDXY", l2PVDXY);
        b.SetBranch("e3PVDXY", l3PVDXY);
        b.SetBranch("e1PVDZ", l1PVDZ);
        b.SetBranch("e2PVDZ", l2PVDZ);
        b.SetBranch("e3PVDZ", l3PVDZ);
        b.SetBranch("Mass", Mass);
        b.SetBranch("e1_e2_Mass", ZMass);
        b.SetBranch("e1Pt", l1Pt);
        b.SetBranch("e2Pt", l2Pt);
        b.SetBranch("e3Pt", l3Pt);
        b.SetBranch("e1Eta", l1Eta);
        b.SetBranch("e2Eta", l2Eta);
        b.SetBranch("e3Eta", l3Eta);
        b.SetBranch("e1Phi", l1Phi);
        b.SetBranch("e2Phi", l2Phi);
        b.SetBranch("e3Phi", l3Phi);
        b.SetBranch("e1Mass", l1Mass);
        b.SetBranch("e2Mass", l2Mass);
        b.SetBranch("e3Mass", l3Mass);
        b.SetBranch("e3MtToMET", l3MtToMET);
        if (isMC_) {
            b.SetBranch("e1GenPt", l1GenPt);
            b.SetBranch("e2GenPt", l2GenPt);
            b.SetBranch("e3GenPt", l3GenPt);
        }
    }
    else if (channel_ == eem) {
        b.SetBranch("e1IsCBVIDTight", l1IsTight);
        b.SetBranch("e2IsCBVIDTight", l2IsTight);
        b.SetBranch("mIsWZTight", l3IsTight);
        b.SetBranch("mPVDXY", l1PVDXY);
        b.SetBranch("e2PVDXY", l2PVDXY);
        b.SetBranch("e3PVDXY", l3PVDXY);
        b.SetBranch("mPVDZ", l1PVDZ);
        b.SetBranch("e2PVDZ", l2PVDZ);
        b.SetBranch("e3PVDZ", l3PVDZ);
        b.SetBranch("e1_e2_Mass", ZMass);
        b.SetBranch("e1Pt", l1Pt);
        b.SetBranch("e2Pt", l2Pt);
        b.SetBranch("mPt", l3Pt);
        b.SetBranch("e1Eta", l1Eta);
        b.SetBranch("e2Eta", l2Eta);
        b.SetBranch("mEta", l3Eta);
        b.SetBranch("e1Phi", l1Phi);
        b.SetBranch("e2Phi", l2Phi);
        b.SetBranch("mPhi", l3Phi);
        b.SetBranch("e1Mass", l1Mass);
        b.SetBranch("e2Mass", l2Mass);
        b.SetBranch("mMass", l3Mass);
        b.SetBranch("mMtToMET", l3MtToMET);
        if (isMC_) {
            b.SetBranch("mGenPt", l3GenPt);
            b.SetBranch("e1GenPt", l1GenPt);
            b.SetBranch("e2GenPt", l2GenPt);
        }
    }
    else if (channel_ == emm) {
        b.SetBranch("eIsCBVIDTight", l3IsTight);
        b.SetBranch("m1IsWZTight", l1IsTight);
        b.SetBranch("m2IsWZTight", l2IsTight);
        b.SetBranch("m1PVDXY", l1PVDXY);
        b.SetBranch("m2PVDXY", l2PVDXY);
        b.SetBranch("ePVDXY", l3PVDXY);
        b.SetBranch("m1PVDZ", l1PVDZ);
        b.SetBranch("m2PVDZ", l2PVDZ);
        b.SetBranch("ePVDZ", l3PVDZ);
        b.SetBranch("m1_m2_Mass", ZMass);
        b.SetBranch("m1Pt", l1Pt);
        b.SetBranch("m2Pt", l2Pt);
        b.SetBranch("ePt", l3Pt);
        b.SetBranch("m1Eta", l1Eta);
        b.SetBranch("m2Eta", l2Eta);
        b.SetBranch("eEta", l3Eta);
        b.SetBranch("m1Phi", l1Phi);
        b.SetBranch("m2Phi", l2Phi);
        b.SetBranch("ePhi", l3Phi);
        b.SetBranch("m1Mass", l1Mass);
        b.SetBranch("m2Mass", l2Mass);
        b.SetBranch("eMass", l3Mass);
        b.SetBranch("eMtToMET", l3MtToMET);
        if (isMC_) {
            b.SetBranch("eGenPt", l3GenPt);
            b.SetBranch("m1GenPt", l1GenPt);
            b.SetBranch("m2GenPt", l2GenPt);
        }
    }
    else if (channel_ == mmm) {
        b.SetBranch("m1IsWZTight", l1IsTight);
        b.SetBranch("m2IsWZTight", l2IsTight);
        b.SetBranch("m3IsWZTight", l3IsTight);
        b.SetBranch("m1PVDXY", l1PVDXY);
        b.SetBranch("m2PVDXY", l2PVDXY);
        b.SetBranch("m3PVDXY", l3PVDXY);
        b.SetBranch("m1PVDZ", l1PVDZ);
        b.SetBranch("m2PVDZ", l2PVDZ);
        b.SetBranch("m3PVDZ", l3PVDZ);
        b.SetBranch("m1_m2_Mass", ZMass);
        b.SetBranch("m1Pt", l1Pt);
        b.SetBranch("m2Pt", l2Pt);
        b.SetBranch("m3Pt", l3Pt);
        b.SetBranch("m1Eta", l1Eta);
        b.SetBranch("m2Eta", l2Eta);
        b.SetBranch("m3Eta", l3Eta);
        b.SetBranch("m1Phi", l1Phi);
        b.SetBranch("m2Phi", l2Phi);
        b.SetBranch("m3Phi", l3Phi);
        b.SetBranch("m1Mass", l1Mass);
        b.SetBranch("m2Mass", l2Mass);
        b.SetBranch("m3Mass", l3Mass);
        b.SetBranch("m3MtToMET", l3MtToMET);
        if (isMC_) {
            b.SetBranch("m1GenPt", l1GenPt);
            b.SetBranch("m2GenPt", l2GenPt);
            b.SetBranch("m3GenPt", l3GenPt);
        }
    }

    b.SetBranch("type1_pfMETEt", MET);
    b.SetBranch("type1_pfMETPhi", type1_pfMETPhi);
    b.SetBranch("nCBVIDTightElec", nCBVIDTightElec);
    b.SetBranch("nCBVIDHLTSafeElec", nCBVIDHLTSafeElec);
    b.SetBranch("nCBVIDVetoElec", nCBVIDVetoElec);
    b.SetBranch("nWZTightMuon", nWZTightMuon);
    b.SetBranch("nWZMediumMuon", nWZMediumMuon);
    b.SetBranch("Flag_BadChargedCandidateFilterPass", Flag_BadChargedCandidateFilterPass);
    b.SetBranch("Flag_BadPFMuonFilterPass", Flag_BadPFMuonFilterPass);
    b.SetBranch("Flag_HBHENoiseFilterPass", Flag_HBHENoiseFilterPass);
    b.SetBranch("Flag_HBHENoiseIsoFilterPass", Flag_HBHENoiseIsoFilterPass);
    b.SetBranch("Flag_EcalDeadCellTriggerPrimitiveFilterPass", Flag_EcalDeadCellTriggerPrimitiveFilterPass);
    b.SetBranch("Flag_goodVerticesPass", Flag_goodVerticesPass);
    b.SetBranch("Flag_eeBadScFilterPass", Flag_eeBadScFilterPass);
    b.SetBranch("Flag_globalTightHalo2016FilterPass", Flag_globalTightHalo2016FilterPass);
}

void WZSelectorBase::SetBranchesNanoAOD() {
    b.CleanUp();
    b.SetBranch("nElectron", nElectron);
    b.SetBranch("nMuon", nMuon);
    b.SetBranch("Electron_pt", Electron_pt);
    b.SetBranch("Electron_eta", Electron_eta);
    b.SetBranch("Electron_phi", Electron_phi);
    b.SetBranch("Electron_mass", Electron_mass);
    b.SetBranch("Muon_pt", Muon_pt);
    b.SetBranch("Muon_eta", Muon_eta);
    b.SetBranch("Muon_phi", Muon_phi);
    b.SetBranch("Muon_mass", Muon_mass);
    b.SetBranch("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
    //b.SetBranch("Electron_cutBased", Electron_cutBased);
    b.SetBranch("Electron_cutBased_Fall17_V1", Electron_cutBased);
    b.SetBranch("Electron_dxy", Electron_dxy);
    b.SetBranch("Electron_dz", Electron_dz);
    b.SetBranch("Muon_tightId", Muon_tightId);
    b.SetBranch("Muon_tkIsoId", Muon_tkIsoId);
    b.SetBranch("Muon_dxy", Muon_dxy);
    b.SetBranch("Muon_dz", Muon_dz);
    b.SetBranch("MET_pt", MET);
    b.SetBranch("MET_phi", type1_pfMETPhi);
    b.SetBranch("Electron_charge", Electron_charge);
    b.SetBranch("Muon_charge", Muon_charge);
    if (isMC_) {
        //b.SetBranch("e1GenPt", l1GenPt);
        //b.SetBranch("e2GenPt", l2GenPt);
        //b.SetBranch("e3GenPt", l3GenPt);
        b.SetBranch("genWeight", genWeight);
        b.SetBranch("Pileup_nPU", numPU);
    }
    //else {
    //    b.SetBranch("Flag_duplicateMuonsPass", Flag_duplicateMuonsPass);
    //    b.SetBranch("Flag_badMuonsPass", Flag_badMuonsPass);
    //}
    //b.SetBranch("Flag_BadChargedCandidateFilterPass", Flag_BadChargedCandidateFilterPass);
    //b.SetBranch("Flag_BadPFMuonFilterPass", Flag_BadPFMuonFilterPass);
    //b.SetBranch("Flag_HBHENoiseFilterPass", Flag_HBHENoiseFilterPass);
    //b.SetBranch("Flag_HBHENoiseIsoFilterPass", Flag_HBHENoiseIsoFilterPass);
    //b.SetBranch("Flag_EcalDeadCellTriggerPrimitiveFilterPass", Flag_EcalDeadCellTriggerPrimitiveFilterPass);
    //b.SetBranch("Flag_goodVerticesPass", Flag_goodVerticesPass);
    //b.SetBranch("Flag_eeBadScFilterPass", Flag_eeBadScFilterPass);
    //b.SetBranch("Flag_globalTightHalo2016FilterPass", Flag_globalTightHalo2016FilterPass);
}

void WZSelectorBase::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    b.SetEntry(entry);
    weight = 1;

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

    l1Pt = 0;
    l2Pt = 0;
    l3Pt = 0;
    l1Eta = 0;
    l2Eta = 0;
    l3Eta = 0;
    l1Phi = 0;
    l2Phi = 0;
    l3Phi = 0;
    l1Mass = 0;
    l2Mass = 0;
    l3Mass = 0;

    SetChannelAndIndicesNano();
    if (channel_ != eee && channel_ != eem && channel_ != emm && channel_ != mmm) {
        passesLeptonVeto = false;
        return;
    }
    else
        passesLeptonVeto = true;

    SetGoodLeptonsFromNano();
    SetLeptonVarsNano();
    SetMasses();

    if (isMC_) {
        ApplyScaleFactors();
    }
    //else {
    //}

}

void WZSelectorBase::SetLeptonVarsNano() {
    auto& l1 = leptons.at(0);
    auto& l2 = leptons.at(1);
    auto& l3 = leptons.at(2);

    l1Pt = l1.pt();
    l1Eta = l1.eta();
    l1Phi = l1.phi();
    l1Mass= l1.mass();
    l2Pt = l2.pt();
    l2Eta = l2.eta();
    l2Phi = l2.phi();
    l2Mass= l2.mass();
    l3Pt = l3.pt();
    l3Eta = l3.eta();
    l3Phi = l3.phi();
    l3Mass= l3.mass();

    if (channel_ == eee) {
        size_t ie = looseElecIndices.at(0);
        l1IsTight = (Electron_cutBased[ie] == 4);
        l1PVDXY = Electron_dxy[ie];
        l1PVDZ = Electron_dz[ie];
        ie = looseElecIndices.at(1);
        l2IsTight = (Electron_cutBased[ie] == 4);
        l2PVDXY = Electron_dxy[ie];
        l2PVDZ = Electron_dz[ie];
        ie = looseElecIndices.at(2);
        l3IsTight = (Electron_cutBased[ie] == 4);
        l3PVDXY = Electron_dxy[ie];
        l3PVDZ = Electron_dz[ie];
    }
    else if (channel_ == eem) {
        size_t ie = looseElecIndices.at(0);
        l1IsTight = (Electron_cutBased[ie] == 4);
        l1PVDXY = Electron_dxy[ie];
        l1PVDZ = Electron_dz[ie];
        ie = looseElecIndices.at(1);
        l2IsTight = (Electron_cutBased[ie] == 4);
        l2PVDXY = Electron_dxy[ie];
        l2PVDZ = Electron_dz[ie];
        size_t imu = looseMuonIndices.at(0);
        l3IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
    }
    else if (channel_ == emm) {
        size_t imu = looseMuonIndices.at(0);
        l1IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
        l1PVDXY = Muon_dxy[imu];
        l1PVDZ = Muon_dz[imu];
        imu = looseMuonIndices.at(1);
        l2IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
        l2PVDXY = Muon_dxy[imu];
        l2PVDZ = Muon_dz[imu];
        size_t ie = looseElecIndices.at(0);
        l3IsTight = (Electron_cutBased[ie] == 4);
        l3PVDXY = Electron_dxy[ie];
        l3PVDZ = Electron_dz[ie];
     }
     else if (channel_ == mmm) {
        size_t imu = looseMuonIndices.at(0);
        l1IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
        l1PVDXY = Muon_dxy[imu];
        l1PVDZ = Muon_dz[imu];
        imu = looseMuonIndices.at(1);
        l2IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
        l2PVDXY = Muon_dxy[imu];
        l2PVDZ = Muon_dz[imu];
        imu = looseMuonIndices.at(2);
        l3IsTight = (Muon_tightId[imu] && Muon_pfRelIso04_all[imu] < 0.15);
        l3PVDXY = Muon_dxy[imu];
        l3PVDZ = Muon_dz[imu];
     }
}

// Always ordered: 0 - Zlep1, 1 - Zlep2, 2 - Wlep (lep IDs are clear from channel)
void WZSelectorBase::SetGoodLeptonsFromNano() {
    leptons.clear();
    bool zee = (channel_ == eem || channel_ == eee);
    auto& indices = zee ? looseElecIndices : looseMuonIndices;
    if (!(indices.size() == 2 || indices.size() == 3)) {
        throw std::length_error("Invalid lepton indices");
    }
    for (const auto& i : indices) {
        auto lep = zee ?
            LorentzVector(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]) :
            LorentzVector(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
        leptons.push_back(lep);
    }

    if (channel_ == eem || channel_ == emm) {
        auto& wIndices = (channel_ == eem) ? looseMuonIndices : looseElecIndices;
        if (wIndices.size() != 1)
            throw std::length_error("Invalid W lepton indices");
        size_t wi = wIndices.at(0);
        auto wlep = (channel_ == eem) ?
                LorentzVector(Muon_pt[wi], Muon_eta[wi], Muon_phi[wi], Muon_mass[wi]) :
                LorentzVector(Electron_pt[wi], Electron_eta[wi], Electron_phi[wi], Electron_mass[wi]);
        leptons.push_back(wlep);
    }
}

void WZSelectorBase::SetChannelAndIndicesNano() {
    // cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
    nCBVIDTightElec = 0;
    nCBVIDVetoElec = 0;
    nWZTightMuon = 0;
    nWZMediumMuon = 0;

    //TODO: Embed these variables in the NanoSkims
    looseMuonIndices.clear();
    looseElecIndices.clear();
    for (size_t i = 0; i < nMuon; i++) {
        if (Muon_tightId[i] && Muon_pfRelIso04_all[i] < 0.40) {
            nWZMediumMuon++;
            looseMuonIndices.push_back(i);
            nWZTightMuon += (Muon_pfRelIso04_all[i] < 0.15);
        }
    }

    for (size_t i = 0; i < nElectron; i++) {
        if (Electron_cutBased[i] >= 1) {
            nCBVIDVetoElec++;
            looseElecIndices.push_back(i);
        }
        if (Electron_cutBased[i] == 4)
            nCBVIDTightElec++;
    }

    if (nWZMediumMuon == 0 && nCBVIDVetoElec == 3)
        channelName_ = "eee";
    else if (nWZMediumMuon == 1 && nCBVIDVetoElec == 2)
        channelName_ = "eem";
    else if (nWZMediumMuon == 2 && nCBVIDVetoElec == 1)
        channelName_ = "emm";
    else if (nWZMediumMuon == 3 && nCBVIDVetoElec == 0)
        channelName_ = "mmm";
    else
        channelName_ = "Unknown";
    
    //std::cout << "Channel " << channelName_ << " " << channel_ << " elecIndices " << looseElecIndices.size() 
    //          << " muon indices " << looseMuonIndices.size() << std::endl;

    channel_ = channelMap_[channelName_];
}

void WZSelectorBase::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation){ 
    b.SetEntry(entry);
    weight = 1;

    if (isMC_) {
        ApplyScaleFactors();
    }
    // Veto on tight leptons
    // Make sure tight leptons also pass loose
    // passesLeptonVeto = nWZTightMuon + nCBVIDTightElec <= 3 &&
    //                      (nWZMediumMuon + nCBVIDHLTSafeElec) >= (nWZTightMuon+nCBVIDTightElec);
    // If tight isn't required to include loose
    // passesLeptonVeto = std::abs(nWZMediumMuon + nCBVIDHLTSafeElec - (l1IsLoose +l2IsLoose +l3IsLoose)) < 0.1;

    // Veto on loose leptons
    //passesLeptonVeto = (nWZMediumMuon + nCBVIDHLTSafeElec) == 3;
    passesLeptonVeto = (nWZMediumMuon + nCBVIDVetoElec) == 3;
 
}

void WZSelectorBase::ApplyScaleFactors() {
    weight = genWeight;
    if (channel_ == eee) {
        if (eIdSF_ != nullptr) {
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        if (eGsfSF_ != nullptr) {
            weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
    }
    else if (channel_ == eem) {
        if (eIdSF_ != nullptr) {
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        if (eGsfSF_ != nullptr) {
            weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        if (mIdSF_ != nullptr) {
            weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        if (mIsoSF_ != nullptr) {
            weight *= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
    }
    else if (channel_ == emm) {
        if (eIdSF_ != nullptr) {
            weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        if (eGsfSF_ != nullptr) {
            weight *= eGsfSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        if (mIdSF_ != nullptr) {
            weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        if (mIsoSF_ != nullptr) {
            weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
    }
    else {
        if (mIdSF_ != nullptr) {
            weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
        if (mIsoSF_ != nullptr) {
            weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
    }
    if (pileupSF_ != nullptr) {
        //weight *= pileupSF_->Evaluate1D(numPU);
        weight *= pileupSF_->Evaluate1D(nTruePU);
    }
}

void WZSelectorBase::SetMasses() {
    if (leptons.size() == 0) {
        leptons.push_back(LorentzVector(l1Pt, l1Eta, l1Phi, l1Mass));
        leptons.push_back(LorentzVector(l2Pt, l2Eta, l2Phi, l2Mass));
        leptons.push_back(LorentzVector(l3Pt, l3Eta, l3Phi, l3Mass));
    }

    ZMass = (leptons.at(0)+leptons.at(1)).M();
    Mass = (leptons.at(0)+leptons.at(1)+leptons.at(2)).M();
}

// Meant to be a wrapper for the tight ID just in case it changes
// To be a function of multiple variables
bool WZSelectorBase::zlep1IsTight() {
    bool l1IsE = (channel_ == eem || channel_ == eee);
    if (l1IsE) { 
        // TODO: This should really be the supercluster eta, I think
        bool isEB = std::abs(l1Eta) > 1.479;
        return l1IsTight && (
            std::abs(l1PVDXY) < (isEB ? 0.05 : 0.1) && std::abs(l1PVDZ) < (isEB ? 0.1 : 0.2));
    }
    return (l1IsTight && std::abs(l1PVDXY) < 0.02 && std::abs(l1PVDZ) < 0.1); 
}

bool WZSelectorBase::zlep2IsTight() {
    bool l2IsE = (channel_ == eem || channel_ == eee);
    if (l2IsE) { 
        // TODO: This should really be the supercluster eta, I think
        bool isEB = std::abs(l2Eta) > 1.479;
        return l2IsTight && (
            std::abs(l2PVDXY) < (isEB ? 0.05 : 0.1) && std::abs(l2PVDZ) < (isEB ? 0.1 : 0.2));
    }
    return (l2IsTight && std::abs(l2PVDXY) < 0.02 && std::abs(l2PVDZ) < 0.1); 
}

bool WZSelectorBase::tightZLeptons() {
    return zlep1IsTight() && zlep2IsTight(); 
}

bool WZSelectorBase::lepton3IsTight() {
    bool l3IsE = (channel_ == emm || channel_ == eee);
    if (l3IsE) { 
        // TODO: This should really be the supercluster eta, I think
        bool isEB = std::abs(l3Eta) > 1.479;
        return l3IsTight && (
            std::abs(l3PVDXY) < (isEB ? 0.05 : 0.1) && std::abs(l3PVDZ) < (isEB ? 0.1 : 0.2));
    }
    return (l3IsTight && std::abs(l3PVDXY) < 0.02 && std::abs(l3PVDZ) < 0.1); 
}

bool WZSelectorBase::IsGenMatched3l() {
    //return true;
    return (!isMC_ || isNonpromptMC_ || 
        (isZgamma_ && l1GenPt > 0 && l2GenPt > 0) ||
        (l1GenPt > 0 && l2GenPt > 0 && l3GenPt > 0));
}
