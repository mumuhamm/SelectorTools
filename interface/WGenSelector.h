#ifndef WGenSelector_h
#define WGenSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <exception>
#include <iostream>

// Headers needed by this particular selector
#include <vector>
#include "Analysis/VVAnalysis/interface/ScaleFactor.h"
#include "Analysis/VVAnalysis/interface/NanoGenSelectorBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class WGenSelector : public NanoGenSelectorBase {
public :
    // Derived values
    LorentzVector nu;
    LorentzVector wCand;
    LorentzVector wCandMet;
    float mWlhe;
    float mTtrue;
    float mTmet;
    float cenWeight;

    float MW_GEN_;
    float GAMMAW_GEN_;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    WGenSelector(TTree * /*tree*/ =0) { }
    ~WGenSelector() { }

    ClassDefOverride(WGenSelector,0);

protected:
    virtual void SetComposite() override;
    virtual void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillHistogramsByName(Long64_t entry, std::string& toAppend, std::pair<Systematic, std::string> variation);
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    double breitWignerWeight(double offset);
};

#endif



