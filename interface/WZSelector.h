#ifndef WZSelector_h
#define WZSelector_h

#include <vector>
#include "Analysis/SelectorTools/interface/WZSelectorBase.h"
#include <TH3.h>

class WZSelector : public WZSelectorBase {
public :
    bool doSystematics_ = false;
    bool applyFullSelection_ = true;
    bool isaQGC_ = false;
    bool doaQGC_ = false;

    std::vector<float>* scaleWeights = NULL;
    std::vector<float>* pdfWeights = NULL;
    std::vector<float> lheWeights;
    unsigned int weight_info_;
    float dEtajj;
    float dEtajj_jesUp;
    float dEtajj_jesDown;
    float dEtajj_jerUp;
    float dEtajj_jerDown;
    float dRjj;
    float dRjj_jesUp;
    float dRjj_jesDown;
    float dRjj_jerUp;
    float dRjj_jerDown;
    float zep3l;
    float zep3l_jesUp;
    float zep3l_jesDown;
    float zep3l_jerUp;
    float zep3l_jerDown;
    std::vector<float>* jetEta = NULL;
    std::vector<float>* jetEta_jesUp = NULL;
    std::vector<float>* jetEta_jesDown = NULL;
    std::vector<float>* jetEta_jerUp = NULL;
    std::vector<float>* jetEta_jerDown = NULL;
    std::vector<float>* jetPhi = NULL;
    std::vector<float>* jetPhi_jesUp = NULL;
    std::vector<float>* jetPhi_jesDown = NULL;
    std::vector<float>* jetPhi_jerUp = NULL;
    std::vector<float>* jetPhi_jerDown = NULL;
    std::vector<float>* jetPt = NULL;
    std::vector<float>* jetPt_jesUp = NULL;
    std::vector<float>* jetPt_jesDown = NULL;
    std::vector<float>* jetPt_jerUp = NULL;
    std::vector<float>* jetPt_jerDown = NULL;
    std::vector<float>* jetCSVv2 = NULL;
    UInt_t nvtx;
    Float_t Zlep1_Wlep_Mass;
    Float_t Zlep2_Wlep_Mass;
    Float_t Eta;
    Float_t Pt;
    Float_t ZPt;
    Float_t ZEta;
    Float_t ZPhi;
    Float_t mjj;
    Float_t mjj_jesUp;
    Float_t mjj_jesDown;
    Float_t mjj_jerUp;
    Float_t mjj_jerDown;
    Float_t MtToMET;
    Float_t MtWZ;
    Float_t M3lMET;
    Float_t type1_pfMETEt_UncTool;
    Float_t type1_pfMETEt_jesUp;
    Float_t type1_pfMETEt_jesDown;
    Float_t type1_pfMETEt_jerUp;
    Float_t type1_pfMETEt_jerDown;
    Float_t type1_pfMETEt_unclusteredEnUp;
    Float_t type1_pfMETEt_unclusteredEnDown;
    Float_t l1PtScaleCorrErr;
    Float_t l2PtScaleCorrErr;
    Float_t l3PtScaleCorrErr;
    
    TBranch* b_l3MtToMET;
    TBranch* b_MtToMET;
    TBranch* b_pdfWeights;
    TBranch* b_scaleWeights;
    TBranch* b_jetCSVv2;
    TBranch* b_jetPt;
    TBranch* b_jetEta;
    TBranch* b_jetEta_jesUp;
    TBranch* b_jetEta_jesDown;
    TBranch* b_jetEta_jerUp;
    TBranch* b_jetEta_jerDown;
    TBranch* b_jetPhi;
    TBranch* b_jetPhi_jesUp;
    TBranch* b_jetPhi_jesDown;
    TBranch* b_jetPhi_jerUp;
    TBranch* b_jetPhi_jerDown;
    TBranch* b_jetPt_jesUp;
    TBranch* b_jetPt_jesDown;
    TBranch* b_jetPt_jerUp;
    TBranch* b_jetPt_jerDown;
    TBranch* b_nvtx;
    TBranch* b_Zlep1_Wlep_Mass;
    TBranch* b_Zlep2_Wlep_Mass;
    TBranch* b_Eta;
    TBranch* b_Pt;
    TBranch* b_ZPt;
    TBranch* b_ZPhi;
    TBranch* b_ZEta;
    TBranch* b_mjj;
    TBranch* b_mjj_jesUp;
    TBranch* b_mjj_jesDown;
    TBranch* b_mjj_jerUp;
    TBranch* b_mjj_jerDown;
    TBranch* b_l1PtScaleCorrErr;
    TBranch* b_l2PtScaleCorrErr;
    TBranch* b_l3PtScaleCorrErr;
    TBranch* b_type1_pfMETEt_jesUp;
    TBranch* b_type1_pfMETEt_UncTool;
    TBranch* b_type1_pfMETEt_jesDown;
    TBranch* b_type1_pfMETEt_jerUp;
    TBranch* b_type1_pfMETEt_jerDown;
    TBranch* b_type1_pfMETEt_unclusteredEnUp;
    TBranch* b_type1_pfMETEt_unclusteredEnDown;

    // Readers to access the data (delete the ones you do not need).
    virtual void    Init(TTree *tree) override;
    virtual void    SetupNewDirectory() override;

    ClassDefOverride(WZSelector,0);
protected:
    virtual void    SetBranchesUWVV() override;
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void FillVBSHistograms(float weight, bool noBlind, 
            std::pair<Systematic, std::string> variation);
    bool PassesBaseSelection(Long64_t entry, bool tightLeps, Selection selection);
    bool PassesVBSSelection(bool noBlind);
    bool PassesVBSBackgroundControlSelection();
    bool PassesFullWZSelection(Long64_t entry);
    unsigned int GetLheWeightInfo();
    void ShiftEfficiencies(Systematic variation);
    float GetMuonScaleUncertainty(float muEta);
};

#endif

