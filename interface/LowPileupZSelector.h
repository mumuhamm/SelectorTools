#ifndef LowPileupZSelector_h
#define LowPileupZSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/SelectorBase.h"
#include "Analysis/VVAnalysis/interface/BranchManager.h"
#include "Analysis/VVAnalysis/interface/LowPileupSelector.h"

#include "Analysis/VVAnalysis/interface/helpers.h"
#include <TSelector.h>

class LowPileupZSelector : public LowPileupSelector {
public :
    // Derived values
    TLorentzVector zCand;

    TLorentzVector* lep1 = nullptr;
    TLorentzVector* lep2 = nullptr;
    Float_t genWeight;
    Float_t PUWeight;

    UInt_t category;

    TBranch* category_b;
    TBranch* lep1_b;
    TBranch* lep2_b;

    enum { eMuMu2HLT=1, eMuMu1HLT, eMuMuNoSel, eMuSta, eMuTrk };
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    LowPileupZSelector(TTree * /*tree*/ =0) { }
    ~LowPileupZSelector() { }

    ClassDefOverride(LowPileupZSelector,0);

protected:
    virtual void    SetBranchesBacon() override;
    void LoadBranchesBacon(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void SetComposite() override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif



