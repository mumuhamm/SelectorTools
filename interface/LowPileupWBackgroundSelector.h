#ifndef LowPileupWBackgroundSelector_h
#define LowPileupWBackgroundSelector_h

#include "Analysis/SelectorTools/interface/LowPileupWSelector.h"

class LowPileupWBackgroundSelector : public LowPileupWSelector {
public :
    LowPileupWBackgroundSelector(TTree * /*tree*/ =0) { }
    ~LowPileupWBackgroundSelector() { }
    virtual void    Init(TTree *tree) override;

    ClassDefOverride(LowPileupWBackgroundSelector,0);

protected:
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void    SetupNewDirectory() override;
};

#endif





