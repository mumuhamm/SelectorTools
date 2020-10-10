#ifndef LowPileupWSelector_h
#define LowPileupWSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <exception>
#include <iostream>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/SelectorTools/interface/ScaleFactor.h"
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include "Analysis/SelectorTools/interface/BranchManager.h"
#include "Analysis/SelectorTools/interface/LowPileupSelector.h"

#include "Analysis/SelectorTools/interface/helpers.h"

class LowPileupWSelector : public LowPileupSelector {
public :
    // Derived values
    TLorentzVector wCand;
    TLorentzVector* met;

    // Mixing quite a lot of things here, but it makes it easier
    TTreeReader     fReader;
    TTreeReaderArray<double> evtWeight = {fReader, "evtWeight"};
    TTreeReaderArray<double> metVector = {fReader, "metVars"};
    TTreeReaderArray<double> metPhiVector = {fReader, "metVarsPhi"};
    TTreeReaderValue<float> uncorrMet = {fReader, "met"};
    TTreeReaderValue<float> uncorrMetPhi = {fReader, "metPhi"};
    TTreeReaderValue<TLorentzVector> lep = {fReader, "lep"};
    TTreeReaderValue<double> mtW = {fReader, "mtCorr"};
    TTreeReaderValue<int> charge = {fReader, "q"};
    TTreeReaderValue<double> lepRelIso = {fReader, "relIso"};
    float mtWcalc;

    std::unordered_map<Systematic, size_t> systematicWeightMap_;
    std::unordered_map<Systematic, size_t> metCorrWeightMap_;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    LowPileupWSelector(TTree * /*tree*/ =0) { }
    ~LowPileupWSelector() { }

    ClassDefOverride(LowPileupWSelector,0);

protected:
    virtual void    SetBranchesBacon() override;
    void LoadBranchesBacon(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void SetComposite() override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif




