#include "Analysis/SelectorTools/interface/LowPileupZSelector.h"
#include <TStyle.h>
#include <regex>

void LowPileupZSelector::Init(TTree *tree)
{
    allChannels_ = {{mm, "mm"}};
    hists1D_ = {"CutFlow", "mZ", "yZ", "ptZ", "ptl1", "etal1", "ptl2", "etal2", "pfMet",};

    LowPileupSelector::Init(tree);
}

void LowPileupZSelector::SetBranchesBacon() {
    channel_ = mm;
    channelName_ = "mm";
    lep1 = nullptr;
    lep2 = nullptr;
    fChain->SetBranchAddress("category", &category, &category_b);
    //fChain->SetBranchAddress("met", &pfMet, &b_pfMet);
    //fChain->SetBranchAddress("metPhi", &pfMetPhi, &b_pfMetPhi);
    fChain->SetBranchAddress("lep1", &lep1, &lep1_b);
    fChain->SetBranchAddress("lep2", &lep2, &lep2_b);
    if (isMC_) {
        //b.SetBranch("genWeight", genWeight);
        //b.SetBranch("PUWeight", PUWeight);
    }
    LowPileupSelector::SetBranchesBacon();
}

void LowPileupZSelector::LoadBranchesBacon(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    lep1_b->GetEntry(entry);
    lep2_b->GetEntry(entry);
    category_b->GetEntry(entry);
    LowPileupSelector::LoadBranchesBacon(entry, variation);
}

void LowPileupZSelector::SetComposite() {
    zCand = *lep1+*lep2;
}

void LowPileupZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    if (zCand.M() < 60 || zCand.M() > 120)
        return;
    if (category != eMuMu1HLT && category != eMuMu2HLT)
        return;
    //SafeHistFill(histMap1D_, getHistName("mZ", variation.second), zCand.M(), weight);
    //SafeHistFill(histMap1D_, getHistName("ptZ", variation.second), zCand.Pt(), weight);
    //SafeHistFill(histMap1D_, getHistName("yZ", variation.second), zCand.Rapidity(), weight);
    //SafeHistFill(histMap1D_, getHistName("ptl1", variation.second), lep1->Pt(), weight);
    //SafeHistFill(histMap1D_, getHistName("ptl2", variation.second), lep2->Pt(), weight);
    //SafeHistFill(histMap1D_, getHistName("etal1", variation.second), lep1->Eta(), weight);
    //SafeHistFill(histMap1D_, getHistName("etal2", variation.second), lep2->Eta(), weight);
    //SafeHistFill(histMap1D_, getHistName("pfMet", variation.second), pfMet, weight);
}
