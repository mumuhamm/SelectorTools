#include "Analysis/VVAnalysis/interface/LowPileupWSelector.h"
#include <TStyle.h>
#include <regex>

void LowPileupWSelector::Init(TTree *tree)
{
    doSystematics_ = false;
    //doSystematics_ = true;
    allChannels_ = {{ep, "ep"}, {en, "en"}, {mp, "mp"}, {mn, "mn"}};

    hists1D_ = {"CutFlow", "mW", "mtW", "mtWUncorr", "yW", "ptW", "phiW", "ptl", "etal", "pfMet", "pfMetPhi",};
    systHists_ = {"ptW", "ptl", "mtW", "pfMet"};
    systematics_ = {
        {muonEfficiencyMCSubtractUp, "CMS_eff_MCsubt_mUp"},
        {muonEfficiencyMCSubtractDown, "CMS_eff_MCsubt_mDown"},
        {modelingFsrUp, "CMS_modeling_fsrUp"},
        {modelingFsrDown, "CMS_modeling_fsrDown"},
        {muonEfficiencyBackgroundUp, "CMS_eff_background_mUp"},
        {muonEfficiencyBackgroundDown, "CMS_eff_background_mDown"},
        {muonEfficiencyTagPtUp, "CMS_eff_tagPt_mUp"},
        {muonEfficiencyTagPtDown, "CMS_eff_tagPt_mDown"},
        {muonEfficiencyStatUp, "CMS_eff_stat_mUp"},
        {muonEfficiencyStatDown, "CMS_eff_stat_mDown"},
        {recoilCorrectionEtaShapeUp, "CMS_recoilCorrection_EtaShapeUp"},
        {recoilCorrectionRUShapesUp, "CMS_recoilCorrection_RUShapeUp"}, 
        {recoilCorrectionKeysShapeUp, "CMS_recoilCorrection_KeysShapeUp"},
        {recoilCorrectionStat0Up, "CMS_recoilCorrection_StatBin0Up"},    
        {recoilCorrectionStat1Up, "CMS_recoilCorrection_StatBin1Up"},    
        {recoilCorrectionStat2Up, "CMS_recoilCorrection_StatBin2Up"},        
        {recoilCorrectionStat3Up, "CMS_recoilCorrection_StatBin3Up"},        
        {recoilCorrectionStat4Up, "CMS_recoilCorrection_StatBin4Up"},        
        {recoilCorrectionStat5Up, "CMS_recoilCorrection_StatBin5Up"},        
        {recoilCorrectionStat6Up, "CMS_recoilCorrection_StatBin6Up"},        
        {recoilCorrectionStat7Up, "CMS_recoilCorrection_StatBin7Up"},        
        {recoilCorrectionStat8Up, "CMS_recoilCorrection_StatBin8Up"},        
        {recoilCorrectionStat9Up, "CMS_recoilCorrection_StatBin9Up"},        
        {recoilCorrectionEtaShapeDown, "CMS_recoilCorrection_EtaShapeDown"},
        {recoilCorrectionRUShapesDown, "CMS_recoilCorrection_RUShapeDown"}, 
        {recoilCorrectionKeysShapeDown, "CMS_recoilCorrection_KeysShapeDown"},
        {recoilCorrectionStat0Down, "CMS_recoilCorrection_StatBin0Down"},    
        {recoilCorrectionStat1Down, "CMS_recoilCorrection_StatBin1Down"},    
        {recoilCorrectionStat2Down, "CMS_recoilCorrection_StatBin2Down"},        
        {recoilCorrectionStat3Down, "CMS_recoilCorrection_StatBin3Down"},        
        {recoilCorrectionStat4Down, "CMS_recoilCorrection_StatBin4Down"},        
        {recoilCorrectionStat5Down, "CMS_recoilCorrection_StatBin5Down"},        
        {recoilCorrectionStat6Down, "CMS_recoilCorrection_StatBin6Down"},        
        {recoilCorrectionStat7Down, "CMS_recoilCorrection_StatBin7Down"},        
        {recoilCorrectionStat8Down, "CMS_recoilCorrection_StatBin8Down"},        
        {recoilCorrectionStat9Down, "CMS_recoilCorrection_StatBin9Down"},        
    };
    
    //enum{mc=1,fsr,bkg,tagpt,effstat,pfireu,pfired};
    systematicWeightMap_ = {
        {Central, 0},
        {muonEfficiencyMCSubtractUp, 1},
        {muonEfficiencyMCSubtractDown, 1},
        {modelingFsrUp, 2},
        {modelingFsrDown, 2},
        {muonEfficiencyBackgroundUp, 3},
        {muonEfficiencyBackgroundDown, 3},
        {muonEfficiencyTagPtUp, 4},
        {muonEfficiencyTagPtDown, 4},
        {muonEfficiencyStatUp, 5},
        {muonEfficiencyStatDown, 5},
    };

    metCorrWeightMap_ = {
        {recoilCorrectionEtaShapeUp, 2},
        {recoilCorrectionRUShapesUp, 3}, 
        {recoilCorrectionKeysShapeUp, 4},
        {recoilCorrectionStat0Up, 5},    
        {recoilCorrectionStat1Up, 6},    
        {recoilCorrectionStat2Up, 7},        
        {recoilCorrectionStat3Up, 8},        
        {recoilCorrectionStat4Up, 9},        
        {recoilCorrectionStat5Up, 10},        
        {recoilCorrectionStat6Up, 11},        
        {recoilCorrectionStat7Up, 12},        
        {recoilCorrectionStat8Up, 13},        
        {recoilCorrectionStat9Up, 14},        
        {recoilCorrectionEtaShapeDown, 2},
        {recoilCorrectionRUShapesDown, 3}, 
        {recoilCorrectionKeysShapeDown, 4},
        {recoilCorrectionStat0Down, 5},    
        {recoilCorrectionStat1Down, 6},    
        {recoilCorrectionStat2Down, 7},        
        {recoilCorrectionStat3Down, 8},        
        {recoilCorrectionStat4Down, 9},        
        {recoilCorrectionStat5Down, 10},        
        {recoilCorrectionStat6Down, 11},        
        {recoilCorrectionStat7Down, 12},        
        {recoilCorrectionStat8Down, 13},        
        {recoilCorrectionStat9Down, 14},        
    };

    fReader.SetTree(tree);
    LowPileupSelector::Init(tree);
}

void LowPileupWSelector::SetBranchesBacon() {
    LowPileupSelector::SetBranchesBacon();
}

void LowPileupWSelector::LoadBranchesBacon(Long64_t entry, SystPair variation) { 
    if (variation.first == Central) {
        fReader.SetLocalEntry(entry);
        LowPileupSelector::LoadBranchesBacon(entry, variation);
        if (*charge > 0) {
            // This is probably not safe at all
            if (isE_) {
                channel_ = ep ;
                channelName_ = "ep";
            }
            else {
                channel_ = mp ;
                channelName_ = "mp";
            }
        }
        else {
            if (isE_) {
                channel_ = en ;
                channelName_ = "en";
            }
            else {
                channel_ = mn ;
                channelName_ = "mn";
            }
        }
    }
    int metIndex = (isMC_ && (isW_ || isZ_) && !isNonprompt_);
    pfMet = metVector.At(metIndex);
    pfMetPhi = metPhiVector.At(metIndex);
    if (!isMC_) {
        pfMet = *uncorrMet + 3.14159;
        pfMetPhi = *uncorrMetPhi + 3.14156;
    }

    if (isMC_) {
        float cenwgt = evtWeight[systematicWeightMap_[Central]];
        float wgt = evtWeight[systematicWeightMap_[variation.first]];
        weight = cenwgt;
        if (variation.first == muonEfficiencyMCSubtractUp ||
                    variation.first == muonEfficiencyBackgroundUp ||
                    variation.first == muonEfficiencyTagPtUp ||
                    variation.first == muonEfficiencyStatUp ||
                    variation.first == modelingFsrUp) {
            weight = wgt > cenwgt ? wgt : cenwgt + std::abs(cenwgt-wgt);
        }
        else if (variation.first == muonEfficiencyMCSubtractDown ||
                    variation.first == muonEfficiencyBackgroundDown ||
                    variation.first == muonEfficiencyTagPtDown ||
                    variation.first == muonEfficiencyStatDown ||
                    variation.first == modelingFsrDown) {
            weight = wgt < cenwgt ? wgt : cenwgt - std::abs(cenwgt-wgt);
        }
        else if (variation.first == recoilCorrectionEtaShapeUp ||
                variation.first == recoilCorrectionRUShapesUp ||
                variation.first == recoilCorrectionKeysShapeUp ||
                variation.first == recoilCorrectionStat0Up || 
                variation.first == recoilCorrectionStat1Up || 
                variation.first == recoilCorrectionStat2Up || 
                variation.first == recoilCorrectionStat3Up || 
                variation.first == recoilCorrectionStat4Up || 
                variation.first == recoilCorrectionStat5Up || 
                variation.first == recoilCorrectionStat6Up || 
                variation.first == recoilCorrectionStat7Up || 
                variation.first == recoilCorrectionStat8Up || 
                variation.first == recoilCorrectionStat9Up) {
            double tempMet = metVector.At(metCorrWeightMap_[variation.first]);
            double tempMetPhi = metPhiVector.At(metCorrWeightMap_[variation.first]);
            pfMet = tempMet > pfMet ? tempMet : pfMet + (pfMet - tempMet);
            pfMetPhi = tempMet > pfMet ? tempMetPhi : pfMetPhi + (pfMetPhi - tempMetPhi);
        }
        else if (variation.first == recoilCorrectionEtaShapeDown ||
                variation.first == recoilCorrectionRUShapesDown ||
                variation.first == recoilCorrectionKeysShapeDown ||
                variation.first == recoilCorrectionStat0Down || 
                variation.first == recoilCorrectionStat1Down || 
                variation.first == recoilCorrectionStat2Down || 
                variation.first == recoilCorrectionStat3Down || 
                variation.first == recoilCorrectionStat4Down || 
                variation.first == recoilCorrectionStat5Down || 
                variation.first == recoilCorrectionStat6Down || 
                variation.first == recoilCorrectionStat7Down || 
                variation.first == recoilCorrectionStat8Down || 
                variation.first == recoilCorrectionStat9Down) {
            double tempMet = metVector.At(metCorrWeightMap_[variation.first]);
            double tempMetPhi = metPhiVector.At(metCorrWeightMap_[variation.first]);
            pfMet = tempMet < pfMet ? tempMet : pfMet - (tempMet - pfMet);
            pfMetPhi = tempMet < pfMet ? tempMetPhi : pfMetPhi - (tempMetPhi - pfMetPhi);
        }
    }
}

void LowPileupWSelector::SetComposite() {
    pfMetVec = TLorentzVector();
    pfMetVec.SetPtEtaPhiM(pfMet, 0., pfMetPhi, 0.);
    wCand = *lep + pfMetVec;
    //mtW = std::sqrt(2*lep->Pt()*pfMet*(1-std::cos(lep->Phi()-pfMetVec.Phi())));
}

void LowPileupWSelector::FillHistograms(Long64_t entry, SystPair variation) { 
    if (lep->Pt() < 25 || *mtW < 40)
        return;
    SafeHistFill(histMap1D_, "mW", channel_, variation.first, wCand.M(), weight);
    SafeHistFill(histMap1D_, "mtW", channel_, variation.first, *mtW, weight);
    SafeHistFill(histMap1D_, "mtWUncorr", channel_, variation.first, *mtWuncorr, weight);
    SafeHistFill(histMap1D_, "ptW", channel_, variation.first, wCand.Pt(), weight);
    SafeHistFill(histMap1D_, "yW", channel_, variation.first, wCand.Rapidity(), weight);
    SafeHistFill(histMap1D_, "phiW", channel_, variation.first, wCand.Phi(), weight);
    SafeHistFill(histMap1D_, "ptl", channel_, variation.first, lep->Pt(), weight);
    SafeHistFill(histMap1D_, "etal", channel_, variation.first, lep->Eta(), weight);
    SafeHistFill(histMap1D_, "pfMet", channel_, variation.first, pfMet, weight);
    SafeHistFill(histMap1D_, "pfMetPhi", channel_, variation.first, pfMetPhi, weight);

    if (subprocessHistMaps1D_.empty()) {
        return;
    }

    //std::vector<int> binning = {0, 13, 26, 38, 50, 62, 75, 100};
    std::vector<int> binning = {0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100};
    size_t upperIndex = std::distance(binning.begin(), std::upper_bound(binning.begin(), binning.end(), genVPt));

    std::string binname = name_.substr(0, name_.size()-3);
    binname.append("_GenPtW_");
    if (upperIndex == binning.size())
        binname.append(std::to_string(binning.back()));
    else {
        binname.append(std::to_string(binning.at(upperIndex-1)));
        binname.append("_");
        binname.append(std::to_string(binning.at(upperIndex)));
    }
    binname.append(isE_ ? "__e" : "__m");

    try {
        HistMap1D& subprocessMap = subprocessHistMaps1D_.at(binname);
        SafeHistFill(subprocessMap, "mW", channel_, variation.first, wCand.M(), weight);
        SafeHistFill(subprocessMap, "ptW", channel_, variation.first, wCand.Pt(), weight);
        SafeHistFill(subprocessMap, "yW", channel_, variation.first, wCand.Rapidity(), weight);
        SafeHistFill(subprocessMap, "ptl", channel_, variation.first, lep->Pt(), weight);
        SafeHistFill(subprocessMap, "etal", channel_, variation.first, lep->Eta(), weight);
        SafeHistFill(subprocessMap, "pfMet", channel_, variation.first, pfMet, weight);
    }
    catch (std::exception& e) {
        throw std::range_error("Could not find bin " + binname + " in subprocessMap!");
    }
}
