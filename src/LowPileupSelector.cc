#include "Analysis/SelectorTools/interface/LowPileupSelector.h"
#include <TStyle.h>
#include <regex>

void LowPileupSelector::Init(TTree *tree)
{
    SelectorBase::Init(tree);
    if (name_.find("wlnu") != std::string::npos){
        isW_ = true;
    }
    else if (name_.find("DYm50") != std::string::npos){
        isZ_ = true;
    }
}

void LowPileupSelector::SetBranchesBacon() {
    if (isMC_ && (fChain->GetListOfBranches()->FindObject("genV") != nullptr)) {
        genV = nullptr;
        fChain->SetBranchAddress("genV", &genV, &b_genV);
    }
}

void LowPileupSelector::LoadBranchesBacon(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    weight = 1;

    if (isMC_) {
        weight = scale1fb*1000;
        //weight = genWeight*PUWeight*scale1fb;
        if (isZ_ || isW_) {
            b_genV->GetEntry(entry);
        }
    }
    SetComposite();
}

void LowPileupSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();
    if (name_.find("__m") != std::string::npos || name_.find("__tm") != std::string::npos) {
        allChannels_ = {{mp, "mp"}, {mn, "mn"}};
        isE_ = false;
    }
    else if (name_.find("__e") != std::string::npos || name_.find("__te") != std::string::npos) {
        allChannels_ = {{ep, "ep"}, {en, "en"}};
        isE_ = true;
    }

    InitializeHistogramsFromConfig();
}

