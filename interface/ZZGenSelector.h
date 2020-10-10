#ifndef ZZGenSelector_h
#define ZZGenSelector_h

#include <vector>
#include "Analysis/SelectorTools/interface/SelectorBase.h"
#include <TH3.h>

class ZZGenSelector : public SelectorBase {
public :
    
    float GendPhiZZ; //DeltaPhi between Z1 and Z2
    Float_t GenMass;
    Float_t GenPt;
    Float_t GenEta;
    Float_t Genweight;
    Float_t GenZ1mass;
    Float_t GenZ2mass;
    Float_t GenZ1pt;
    Float_t GenZ2pt;
    Float_t GenZ1Phi;
    Float_t GenZ2Phi;
 
    Float_t Genl1Pt;
    Float_t Genl2Pt;
    Float_t Genl3Pt;
    Float_t Genl4Pt;
    Float_t Genl1Eta;
    Float_t Genl2Eta;
    Float_t Genl3Eta;
    Float_t Genl4Eta;
    Float_t Genl1Phi;
    Float_t Genl2Phi;
    Float_t Genl3Phi;
    Float_t Genl4Phi;
    
    TBranch* b_GenZ1mass;
    TBranch* b_GenZ2mass;
    TBranch* b_GenZ1pt;
    TBranch* b_GenZ2pt;
    TBranch* b_GenZ1Phi;
    TBranch* b_GenZ2Phi;
    TBranch* b_GenMass;
    TBranch* b_GenPt;
    TBranch* b_GenEta;

    TBranch* b_Genl1Pt;
    TBranch* b_Genl2Pt;
    TBranch* b_Genl3Pt;
    TBranch* b_Genl4Pt;
    TBranch* b_Genl1Eta;
    TBranch* b_Genl2Eta;
    TBranch* b_Genl3Eta;
    TBranch* b_Genl4Eta;
    TBranch* b_Genl1Phi;
    TBranch* b_Genl2Phi;
    TBranch* b_Genl3Phi;
    TBranch* b_Genl4Phi; 

    // Readers to access the data (delete the ones you do not need).
    ZZGenSelector(TTree * /*tree*/ =0) { }
    virtual ~ZZGenSelector() { }
    virtual void    Init(TTree *tree) override;
    virtual void    SetupNewDirectory() override;

    ClassDefOverride(ZZGenSelector,0);
protected:
    void LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    void LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    virtual void    SetBranchesNanoAOD() override;
    virtual void    SetBranchesUWVV() override;
    
    void SetVariables(Long64_t entry); 
    bool ZZSelection();
    bool Z4lSelection();
    bool ZSelection();
    void FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) override;
    bool e1e2IsZ1();
};

#endif

