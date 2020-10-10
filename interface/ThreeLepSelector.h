#ifndef ThreeLepSelector_h
#define ThreeLepSelector_h

// #include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

// Headers needed by this particular selector
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include "Analysis/SelectorTools/interface/BranchManager.h"
#include "Analysis/SelectorTools/interface/GoodParticle.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;

enum PID {PID_MUON = 13, PID_ELECTRON = 11, PID_BJET = 5, PID_CJET = 4, PID_JET};

class ThreeLepSelector : public SelectorBase {
public :
   #include "Analysis/SelectorTools/interface/FourTopScales.h"
    
    /*****************************************/
    /* ____  ____   ___  __  __   ___ __  __ */
    /* || )) || \\ // \\ ||\ ||  //   ||  || */
    /* ||=)  ||_// ||=|| ||\\|| ((    ||==|| */
    /* ||_)) || \\ || || || \||  \\__ ||  || */
    /*****************************************/

    ScaleFactor* pileupSF_;
    ScaleFactor* muonSF_;
    ScaleFactor* eIdSF_ ;
    ScaleFactor* eGsfSF_;
    ScaleFactor* mIdSF_;
    ScaleFactor* mIsoSF_;

    // Common variables
    Float_t genWeight;
    Float_t MET;
    Float_t type1_pfMETPhi;

    //NanoAOD variables
    static const unsigned int N_KEEP_MU_E_ = 15;
    static const unsigned int N_KEEP_JET_ = 35;
    static const unsigned int N_KEEP_GEN_ = 300;

    
    UInt_t    nElectron;
    Float_t   Electron_pt[N_KEEP_MU_E_];
    Float_t   Electron_eta[N_KEEP_MU_E_];
    Float_t   Electron_phi[N_KEEP_MU_E_];
    Float_t   Electron_mass[N_KEEP_MU_E_];
    Int_t     Electron_cutBased[N_KEEP_MU_E_];
    Int_t     Electron_charge[N_KEEP_MU_E_];
    Float_t   Electron_MVA[N_KEEP_MU_E_];
    Float_t   Electron_miniPFRelIso_all[N_KEEP_MU_E_];
    Float_t   Electron_dxy[N_KEEP_MU_E_];
    Float_t   Electron_dz[N_KEEP_MU_E_];
    Float_t   Electron_sip3d[N_KEEP_MU_E_];
    Bool_t    Electron_convVeto[N_KEEP_MU_E_];
    UChar_t   Electron_lostHits[N_KEEP_MU_E_];
    Int_t     Electron_tightCharge[N_KEEP_MU_E_];
    Float_t   Electron_sieie[N_KEEP_MU_E_];
    Float_t   Electron_hoe[N_KEEP_MU_E_];
    Float_t   Electron_deltaEtaSC[N_KEEP_MU_E_];
    Float_t   Electron_eInvMinusPInv[N_KEEP_MU_E_];
    Float_t   Electron_dr03EcalRecHitSumEt[N_KEEP_MU_E_];
    Float_t   Electron_dr03HcalDepth1TowerSumEt[N_KEEP_MU_E_];
    Float_t   Electron_dr03TkSumPt[N_KEEP_MU_E_];
    Int_t     Electron_vidBitmap[N_KEEP_MU_E_];
    Float_t   Electron_jetRelIso[N_KEEP_MU_E_];
    Float_t   Electron_eCorr[N_KEEP_MU_E_];
  
    UInt_t    nMuon;
    Float_t   Muon_pt[N_KEEP_MU_E_];
    Float_t   Muon_eta[N_KEEP_MU_E_];
    Float_t   Muon_phi[N_KEEP_MU_E_];
    Float_t   Muon_mass[N_KEEP_MU_E_];
    Int_t     Muon_charge[N_KEEP_MU_E_];
    Bool_t    Muon_tightId[N_KEEP_MU_E_];
    Bool_t    Muon_mediumId[N_KEEP_MU_E_];
    UChar_t   Muon_tkIsoId[N_KEEP_MU_E_];
    Float_t   Muon_pfRelIso04_all[N_KEEP_MU_E_];
    Float_t   Muon_miniPFRelIso_all[N_KEEP_MU_E_];
    Float_t   Muon_dxy[N_KEEP_MU_E_];
    Float_t   Muon_dz[N_KEEP_MU_E_];
    Float_t   Muon_sip3d[N_KEEP_MU_E_];
    Bool_t    Muon_isGlobal[N_KEEP_MU_E_];
    Bool_t    Muon_isTracker[N_KEEP_MU_E_];
    Bool_t    Muon_isPFcand[N_KEEP_MU_E_];
    Int_t     Muon_tightCharge[N_KEEP_MU_E_];
    Float_t   Muon_jetPtRelv2[N_KEEP_MU_E_];
    Float_t   Muon_jetRelIso[N_KEEP_MU_E_];

  
    Int_t     numPU;
    Float_t   Pileup_nTrueInt;

    UInt_t    nJet;
    Float_t   Jet_btagCSVV2[N_KEEP_JET_];
    Float_t   Jet_btagDeepB[N_KEEP_JET_];
    Float_t   Jet_eta[N_KEEP_JET_];
    Float_t   Jet_phi[N_KEEP_JET_];
    Float_t   Jet_pt[N_KEEP_JET_];
    Float_t   Jet_mass[N_KEEP_JET_];
    Float_t   Jet_neHEF[N_KEEP_JET_];
    Float_t   Jet_neEmEF[N_KEEP_JET_];
    Int_t     Jet_nConstituents[N_KEEP_JET_];
    Float_t   Jet_chHEF[N_KEEP_JET_];  
    Float_t   Jet_chEmEF[N_KEEP_JET_];
    Int_t     Jet_jetId[N_KEEP_JET_];
    Int_t     Jet_hadronFlavour[N_KEEP_JET_];
    Float_t   Jet_rawFactor[N_KEEP_JET_];
    Float_t   Jet_L1[N_KEEP_JET_];
    Float_t   Jet_L2L3[N_KEEP_JET_];
    
    Int_t GenPart_pdgId[N_KEEP_GEN_];
    Int_t GenPart_genPartIdxMother[N_KEEP_GEN_];
    Int_t GenPart_status[N_KEEP_GEN_];
    UInt_t nGenPart;

    Bool_t Flag_goodVertices;
    Bool_t Flag_globalSuperTightHalo2016Filter;
    Bool_t Flag_HBHENoiseFilter;
    Bool_t Flag_HBHENoiseIsoFilter;
    Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter;
    Bool_t Flag_BadPFMuonFilter;
    Bool_t Flag_ecalBadCalibFilter;
    
    
    ClassDefOverride(ThreeLepSelector,0);

    /*******************************************************/
    /* __ __  ___  ____  __  ___  ____  __     ____  __    */
    /* || || // \\ || \\ || // \\ || )) ||    ||    (( \   */
    /* \\ // ||=|| ||_// || ||=|| ||=)  ||    ||==   \\    */
    /*  \V/  || || || \\ || || || ||_)) ||__| ||___ \_))   */
    /*******************************************************/
  
    Float_t weight;
    BranchManager b;
    std::vector<GoodPart> goodLeptons;
    std::vector<GoodPart> looseLeptons;
    std::vector<GoodPart> goodJets;
    double HT;
    int nJets, nBJets;
    bool passZVeto;
    BTagCalibration calib;
    BTagCalibrationReader btag_reader; // central sys type
    TH2D* Beff_b;
    TH2D* Beff_j;
    ULong64_t event;
    UInt_t lumi;

    TH2D* h_btag_eff_b;
    TH2D* h_btag_eff_c;
    TH2D* h_btag_eff_udsg;

    
    std::map<int, bool> eventVec;
    std::map<int, std::string> info;

    Bool_t HLT_DoubleMu8_Mass8_PFHT300;
    Bool_t HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
    Bool_t HLT_AK8PFJet450;
    Bool_t HLT_PFJet450;

    Float_t PV_x;
    Float_t PV_y;
    Float_t PV_z;
  
    /************************************************************/
    /* _____ __ __ __  __   ___ ______ __   ___   __  __  __    */
    /* ||    || || ||\ ||  //   | || | ||  // \\  ||\ || (( \   */
    /* ||==  || || ||\\|| ((      ||   || ((   )) ||\\||  \\    */
    /* ||    \\_// || \||  \\__   ||   ||  \\_//  || \|| \_))   */
    /************************************************************/

    void setupMuons();
    void setupElectrons();
    void setupJets();
    void setupChannel();
    void printInfo();
  
    bool isGoodMuon(size_t);
    bool isLooseMuon(size_t);
    bool isGoodJet(size_t);
    bool isGoodBJet(size_t);
    bool isGoodElectron(size_t);
    bool isLooseElectron(size_t);

    bool passMVACut(std::vector<std::vector<double> >, int);
    double mvaInterpolate(double pt, std::vector<double> );
    std::map<int, std::vector<std::vector<double> > >  mvaValues;

    
    bool isLooseMVAElectron(size_t);

    size_t getCloseJetIndex(LorentzVector&, double minDR=10);
    bool doesNotOverlap(size_t);
    bool passFullIso(LorentzVector&, double, double);
    bool doesPassZVeto(GoodPart&, std::vector<GoodPart>&);
    bool passTriggerEmu(size_t);
    double LepRelPt(LorentzVector&, LorentzVector&);
    LorentzVector get4Vector(PID, int);
    bool passFakeableCuts(GoodPart&);
    bool MetFilter();
    float getBtagEffFromFile(double, double, int);
    double getWDecayScaleFactor();
    
    //// General Functions
    int getSRBin() const;
    void clearValues();
    void ApplyScaleFactors();

    // Overloaded or necesary functions
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetupNewDirectory() override;
    // Readers to access the data (delete the ones you do not need).
    virtual void    SetScaleFactors() override;
    virtual void    Init(TTree *tree) override;

    ///ignore
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override {return;}
    virtual void    SetBranchesUWVV() override {return;}
};

#endif

