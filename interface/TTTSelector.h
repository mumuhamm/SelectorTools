#ifndef TTTSelector_h
#define TTTSelector_h

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
#include "Analysis/SelectorTools/interface/BranchManager.h"
#include "Analysis/SelectorTools/interface/GoodParticle.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;

class TTTSelector : public SelectorBase {
 public :
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

  UInt_t nElectron;
  Float_t Electron_pt[N_KEEP_MU_E_];
  Float_t Electron_eta[N_KEEP_MU_E_];
  Float_t Electron_phi[N_KEEP_MU_E_];
  Float_t Electron_mass[N_KEEP_MU_E_];
  Int_t Electron_cutBased[N_KEEP_MU_E_];
  Int_t Electron_charge[N_KEEP_MU_E_];
  Float_t Electron_MVA[N_KEEP_MU_E_];
  Float_t Electron_miniPFRelIso_all[N_KEEP_MU_E_];
  Float_t Electron_dxy[N_KEEP_MU_E_];
  Float_t Electron_dz[N_KEEP_MU_E_];
  Float_t Electron_sip3d[N_KEEP_MU_E_];

  
  UInt_t nMuon;
  Float_t Muon_pt[N_KEEP_MU_E_];
  Float_t Muon_eta[N_KEEP_MU_E_];
  Float_t Muon_phi[N_KEEP_MU_E_];
  Float_t Muon_mass[N_KEEP_MU_E_];
  Int_t Muon_charge[N_KEEP_MU_E_];
  Bool_t Muon_tightId[N_KEEP_MU_E_];
  Bool_t Muon_mediumId[N_KEEP_MU_E_];
  UChar_t Muon_tkIsoId[N_KEEP_MU_E_];
  Float_t Muon_pfRelIso04_all[N_KEEP_MU_E_];
  Float_t Muon_miniPFRelIso_all[N_KEEP_MU_E_];
  Float_t Muon_dxy[N_KEEP_MU_E_];
  Float_t Muon_dz[N_KEEP_MU_E_];
  Float_t Muon_sip3d[N_KEEP_MU_E_];

  Int_t numPU;

  UInt_t nJet;
  Float_t Jet_btagCSVV2[N_KEEP_JET_];
  Float_t Jet_btagDeepB[N_KEEP_JET_];
  Float_t Jet_eta[N_KEEP_JET_];
  Float_t Jet_phi[N_KEEP_JET_];
  Float_t Jet_pt[N_KEEP_JET_];
  Float_t Jet_mass[N_KEEP_JET_];
  Float_t Jet_neHEF[N_KEEP_JET_];
  Float_t Jet_neEmEF[N_KEEP_JET_];
  Int_t Jet_nConstituents[N_KEEP_JET_];
  Float_t Jet_chHEF[N_KEEP_JET_];  
  Float_t Jet_chEmEF[N_KEEP_JET_]; 

  ClassDefOverride(TTTSelector,0);

  /*******************************************************/
  /* __ __  ___  ____  __  ___  ____  __     ____  __    */
  /* || || // \\ || \\ || // \\ || )) ||    ||    (( \   */
  /* \\ // ||=|| ||_// || ||=|| ||=)  ||    ||==   \\    */
  /*  \V/  || || || \\ || || || ||_)) ||__| ||___ \_))   */
  /*******************************************************/
  
  Float_t weight;
  BranchManager b;
  std::vector<GoodPart> goodParts;
  std::vector<int> charge;
  double HT;
  int nTightJet, nBJets;

  /************************************************************/
  /* _____ __ __ __  __   ___ ______ __   ___   __  __  __    */
  /* ||    || || ||\ ||  //   | || | ||  // \\  ||\ || (( \   */
  /* ||==  || || ||\\|| ((      ||   || ((   )) ||\\||  \\    */
  /* ||    \\_// || \||  \\__   ||   ||  \\_//  || \|| \_))   */
  /************************************************************/
  
  /// Captial I for particle definition
  bool IsGoodMuon(size_t);
  bool IsGoodCBElectron(size_t);
  bool IsGoodJet(size_t);
  bool IsGoodBJet(size_t);
  bool IsGoodMVAElectron2016(size_t);
  bool IsGoodMVAElectron2017(size_t);

  // lowercase I for helper function (kinda shitty convention, may change)
  bool isLooseJetId(size_t);
  bool isTightJetId(size_t);
  bool isOverlap(size_t);
  bool passFullIso(LorentzVector&, int, int);

  //// General Functions
  int getSRBin();
  void clearValues();
  void ApplyScaleFactors();

  // Overloaded or necesary functions
  virtual void    SetBranchesNanoAOD() override;
  void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
  void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
  virtual void    SetupNewDirectory() override;
  // Readers to access the data (delete the ones you do not need).
  virtual void    SlaveBegin(TTree *tree) override;
  virtual void    Init(TTree *tree) override;

  ///ignore
  void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
  virtual void    SetBranchesUWVV() override;
};

#endif

