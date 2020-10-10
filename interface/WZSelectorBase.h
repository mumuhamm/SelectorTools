#ifndef WZSelectorBase_h
#define WZSelectorBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include "Analysis/SelectorTools/interface/helpers.h"
#include "Analysis/SelectorTools/interface/BranchManager.h"

class WZSelectorBase : public SelectorBase {
public :
    ScaleFactor* pileupSF_;
    ScaleFactor* muonSF_;
    ScaleFactor* eIdSF_ ;
    ScaleFactor* eGsfSF_;
    ScaleFactor* mIdSF_;
    ScaleFactor* mIsoSF_;

    // Derived variables
    bool isVBS_;
    bool passesLeptonVeto;
    Float_t weight;
    std::vector<LorentzVector> leptons;
    std::vector<LorentzVector> jets;
    std::vector<size_t> looseMuonIndices;
    std::vector<size_t> looseElecIndices;

    BranchManager b;

    // Common variables
    UInt_t nCBVIDVetoElec;
    Float_t genWeight;
    Float_t ZMass;
    Float_t Mass;
    Float_t MET;
    Float_t type1_pfMETPhi;
    Float_t l1Phi;
    Float_t l2Phi;
    Float_t l3Phi;
    Float_t l1Mass;
    Float_t l2Mass;
    Float_t l3Mass;

    //NanoAOD variables
    static const unsigned int N_KEEP_MU_E_ = 20;
    UInt_t nElectron;
    UInt_t nMuon;
    Float_t Electron_pt[N_KEEP_MU_E_];
    Float_t Electron_eta[N_KEEP_MU_E_];
    Float_t Electron_phi[N_KEEP_MU_E_];
    Float_t Electron_mass[N_KEEP_MU_E_];
    Float_t Electron_dxy[N_KEEP_MU_E_];
    Float_t Electron_dz[N_KEEP_MU_E_];
    Int_t Electron_cutBased[N_KEEP_MU_E_];
    Int_t Electron_charge[N_KEEP_MU_E_];
    Float_t Muon_pt[N_KEEP_MU_E_];
    Float_t Muon_eta[N_KEEP_MU_E_];
    Float_t Muon_phi[N_KEEP_MU_E_];
    Float_t Muon_mass[N_KEEP_MU_E_];
    Float_t Muon_dxy[N_KEEP_MU_E_];
    Float_t Muon_dz[N_KEEP_MU_E_];
    Float_t Muon_pfRelIso04_all[N_KEEP_MU_E_];
    Int_t Muon_charge[N_KEEP_MU_E_];
    Bool_t Muon_tightId[N_KEEP_MU_E_];
    Bool_t Muon_mediumId[N_KEEP_MU_E_];
    UChar_t Muon_tkIsoId[N_KEEP_MU_E_];
    Int_t numPU;

    // UWVV variables
    Float_t nTruePU;
    UInt_t nCBVIDTightElec;
    UInt_t nCBVIDHLTSafeElec;
    UInt_t nWZTightMuon;
    UInt_t nWZMediumMuon;
    Bool_t l1IsTight;
    Bool_t l2IsTight;
    Bool_t l3IsTight;
    Float_t l1PVDXY;
    Float_t l2PVDXY;
    Float_t l3PVDXY;
    Float_t l1PVDZ;
    Float_t l2PVDZ;
    Float_t l3PVDZ;
    Float_t l1GenPt;
    Float_t l2GenPt;
    Float_t l3GenPt;
    Bool_t l1IsEB;
    Bool_t l2IsEB;
    Bool_t l3IsEB;
    Float_t m1RelPFIsoDBR04;
    Float_t m2RelPFIsoDBR04;
    Float_t m3RelPFIsoDBR04;
    Float_t l1Eta;
    Float_t l2Eta;
    Float_t l3Eta;
    Float_t l1Pt;
    Float_t l2Pt;
    Float_t l3Pt;
    Float_t l3MtToMET;
    Float_t Flag_BadChargedCandidateFilterPass;
    Float_t Flag_BadPFMuonFilterPass;
    Bool_t Flag_duplicateMuonsPass;
    Bool_t Flag_badMuonsPass;
    Bool_t Flag_HBHENoiseFilterPass;
    Bool_t Flag_HBHENoiseIsoFilterPass;
    Bool_t Flag_EcalDeadCellTriggerPrimitiveFilterPass;
    Bool_t Flag_goodVerticesPass;
    Bool_t Flag_eeBadScFilterPass;
    Bool_t Flag_globalTightHalo2016FilterPass;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    SetScaleFactors() override;
    virtual void    Init(TTree *tree) override;

    ClassDefOverride(WZSelectorBase,0);

protected:
    std::vector<std::string> nonprompt3l_ = {
        "tt-lep", "st-schan", "st-tchan-t", "st-tchan-tbar",
        "st-tw", "st-tbarw", "DYm50", "DYm50-1j",
        "DYm50-2j","DYm50-3j","DYm50-4j", "DYm50__LO",
    };

    bool isNonpromptEstimate_;
    bool isNonpromptMC_;
    bool isZgamma_;
    const float FR_MAX_PT_ = 50;
    const float FR_MAX_ETA_ = 2.5;
    bool zlep1IsTight();
    bool zlep2IsTight();
    bool lepton3IsTight();
    bool tightZLeptons();
    bool IsGenMatched3l();
    virtual std::string GetNameFromFile() override;
    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetBranchesUWVV() override;
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void ApplyScaleFactors();
    void SetLeptonVarsNano();
    void SetMasses();
    void SetGoodLeptonsFromNano();
    void SetChannelAndIndicesNano();
};

#endif

