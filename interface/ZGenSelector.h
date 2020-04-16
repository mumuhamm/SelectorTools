#ifndef ZGenSelector_h
#define ZGenSelector_h

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

class ZGenSelector : public NanoGenSelectorBase {
public :
    // Derived values
    LorentzVector zCand;
    TTreeReaderValue<Int_t> Generator_id1 = {fReader, "Generator_id1"};
    TTreeReaderValue<Int_t> Generator_id2 = {fReader, "Generator_id2"};
    float cenWeight;
    
    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    ZGenSelector(TTree * /*tree*/ =0) { }
    ~ZGenSelector() { }

    ClassDefOverride(ZGenSelector,0);

protected:
    virtual void SetComposite() override;
    virtual void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
};

#endif




