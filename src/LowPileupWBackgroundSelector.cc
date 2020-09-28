#include "Analysis/SelectorTools/interface/LowPileupWBackgroundSelector.h"
#include <TStyle.h>
#include <regex>

void LowPileupWBackgroundSelector::Init(TTree *tree) {
    isNonprompt_ = true;
    LowPileupWSelector::Init(tree);
}

void LowPileupWBackgroundSelector::FillHistograms(Long64_t entry, SystPair variation) { 
    if ((!isE_ && *lepRelIso > 0.25 && *lepRelIso < 0.35) ||
            (isE_ && *lepRelIso > 0.3 && *lepRelIso < 0.45))
        return;
    //weight *= (isE_ ? 2. : 0.13);
    LowPileupWSelector::FillHistograms(entry, variation);
}

void LowPileupWBackgroundSelector::SetupNewDirectory() {
    LowPileupWSelector::SetupNewDirectory();

    for (const auto && obj : *currentHistDir_) {
        std::string name = obj->GetName();
        TNamed* named = dynamic_cast<TNamed*>(obj);
        named->SetName(name.insert(name.length()-2, "Fakes_").c_str());
    }
    
}
