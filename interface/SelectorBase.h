#ifndef SelectorBase_h
#define SelectorBase_h

#include <TROOT.h>
#include <string.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TEfficiency.h>
#include <exception>
#include <iostream>
#include <set>

// Headers needed by this particular selector
#include <vector>
#include <unordered_map>
#include "Analysis/SelectorTools/interface/ScaleFactor.h"

#define PAIR(NAME_) {#NAME_, NAME_}

enum Channel {
    e,           m,         
    ep,          en,        mp,     mn,
    ee,          em,        mm,     
    eee,         eem,       emm,    mmm,
    eeee,        eemm,      mmee,   mmmm,
    Inclusive,   Unknown,   lll, all,
};
  
enum Systematic {
    Central=0,
    jetEnergyScaleUp,          jetEnergyScaleDown,
    jetEnergyResolutionUp,     jetEnergyResolutionDown,
    metUnclusteredEnergyUp,    metUnclusteredEnergyDown,
    muonEfficiencyUp,          muonEfficiencyDown,
    muonScaleUp,               muonScaleDown,
    electronEfficiencyUp,      electronEfficiencyDown,
    electronScaleUp,           electronScaleDown,
    pileupUp,                  pileupDown,
    muonEfficiencyMCSubtractUp, muonEfficiencyMCSubtractDown, 
    modelingFsrUp,             modelingFsrDown, 
    muonEfficiencyBackgroundUp, muonEfficiencyBackgroundDown, 
    muonEfficiencyTagPtUp,     muonEfficiencyTagPtDown, 
    muonEfficiencyStatUp,      muonEfficiencyStatDown, 
    recoilCorrectionEtaShapeUp,   recoilCorrectionEtaShapeDown,  
    recoilCorrectionRUShapesUp,   recoilCorrectionRUShapesDown,  
    recoilCorrectionKeysShapeUp,  recoilCorrectionKeysShapeDown, 
    recoilCorrectionStat0Up,      recoilCorrectionStat0Down,     
    recoilCorrectionStat1Up,      recoilCorrectionStat1Down,
    recoilCorrectionStat2Up,      recoilCorrectionStat2Down,     
    recoilCorrectionStat3Up,      recoilCorrectionStat3Down,     
    recoilCorrectionStat4Up,      recoilCorrectionStat4Down,     
    recoilCorrectionStat5Up,      recoilCorrectionStat5Down,     
    recoilCorrectionStat6Up,      recoilCorrectionStat6Down,     
    recoilCorrectionStat7Up,      recoilCorrectionStat7Down,     
    recoilCorrectionStat8Up,      recoilCorrectionStat8Down,     
    recoilCorrectionStat9Up,      recoilCorrectionStat9Down,     
    BareLeptons, PreFSRLeptons, BornParticles, LHEParticles,
    mWShift100MeVUp, mWShift50MeVUp, mWShift25MeVUp, mWShift20MeVUp, mWShift10MeVUp, 
    mWShift100MeVDown, mWShift50MeVDown, mWShift25MeVDown, mWShift20MeVDown, mWShift10MeVDown, 
    mZShift100MeVUp, mZShift50MeVUp, mZShift25MeVUp, mZShift20MeVUp, mZShift10MeVUp, 
    mZShift100MeVDown, mZShift50MeVDown, mZShift25MeVDown, mZShift20MeVDown, mZShift10MeVDown, 
    ptV0to3, ptV3to5, ptV5to7, ptV7to9, ptV9to12, ptV12to15, ptV15to20, ptV20to27, ptV27to40, ptV40toInf,
    ptV0to3_lhe, ptV3to5_lhe, ptV5to7_lhe, ptV7to9_lhe, ptV9to12_lhe, ptV12to15_lhe, 
    ptV15to20_lhe, ptV20to27_lhe, ptV27to40_lhe, ptV40toInf_lhe,
}; 

    
struct HistLabel {
    std::string name;
    Channel channel;
    Systematic variation;

    bool operator==(const HistLabel& h) const {
        return (name == h.name &&
            channel == h.channel &&
            variation == h.variation);
    };
};

namespace std
{
    template <>
    struct hash<HistLabel>
    {
        size_t operator()(const HistLabel& h) const
        {
            return (std::hash<std::string>()(h.name) ^ std::hash<size_t>()(h.channel) ^
                std::hash<size_t>()(h.variation));

        }
    };
}

class SelectorBase : public TSelector {
 public :
    std::map<std::string, ScaleFactor*> scaleFactors;
    TEfficiency* prefireEff_;
    
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    /*********************************/
    /*  _____ _   _ _   _ __  __      */
    /* | ____| \ | | | | |  \/  |___  */
    /* |  _| |  \| | | | | |\/| / __| */
    /* | |___| |\  | |_| | |  | \__ \ */
    /* |_____|_| \_|\___/|_|  |_|___/ */
    /*********************************/
                              
    enum NtupleType {
        NanoAOD,
        UWVV,
        Bacon,
    };

    enum Selection {
        Default,                      None,
        tightleptons,                 ZZGenFiducial,
        Wselection,                   Zselection,
        WselectionAntiIso,            ZselectionTight,       
        Wselection_Full,              FakeRateSelectionLoose,
        FakeRateSelectionTight,       VBSselection_Loose,
        VBSselection_NoZeppenfeld,    VBSselection_Tight,
        VBSselection_Loose_Full,      VBSselection_NoZeppenfeld_Full,
        VBSselection_Tight_Full,      VBSBackgroundControl,
        VBSBackgroundControlATLAS,    VBSBackgroundControl_Full,
        VBSBackgroundControlLoose,    VBSBackgroundControlLoose_Full,
        Inclusive2Jet,                Inclusive2Jet_Full,
        TightWithLooseVeto,           FourTopPlots,
        FourTopCutBasedEl,            FourTopMVAEl,
        BEfficiency,                  test,
    };

    enum Year {
        yrdefault,      yr2016,      yr2017,      yr2018
    };

    typedef std::pair<Channel, std::string> ChannelPair;

    typedef std::unordered_map<HistLabel, TH1D*> HistMap1D;
    typedef std::unordered_map<HistLabel, TH2D*> HistMap2D;
    typedef std::unordered_map<HistLabel, TH3D*> HistMap3D;
    typedef std::pair<Systematic, std::string> SystPair;
    typedef std::map<Systematic, std::string> SystMap;

    /****************************/
    /*  __  __                  */
    /* |  \/  | __ _ _ __  ___  */
    /* | |\/| |/ _` | '_ \/ __| */
    /* | |  | | (_| | |_) \__ \ */
    /* |_|  |_|\__,_| .__/|___/ */
    /*              |_|         */
    /****************************/

    std::map<std::string, Selection> selectionMap_ = {
        {"Default", Default},
        {"None", None},
        {"tightleptons", tightleptons},
        {"ZZGenFiducial", ZZGenFiducial},
        {"Wselection", Wselection},
        {"WselectionAntiIso", WselectionAntiIso},
        {"Zselection", Zselection},
        {"ZselectionTight", ZselectionTight},
        {"Wselection_Full", Wselection_Full},
        {"FakeRateSelectionLoose", FakeRateSelectionLoose},
        {"FakeRateSelectionTight", FakeRateSelectionTight},
        {"VBSselection_Loose", VBSselection_Loose},
        {"VBSselection_NoZeppenfeld", VBSselection_NoZeppenfeld},
        {"VBSselection_Tight", VBSselection_Tight},
        {"VBSselection_Loose_Full", VBSselection_Loose_Full},
        {"VBSselection_NoZeppenfeld_Full", VBSselection_NoZeppenfeld_Full},
        {"VBSselection_Tight_Full", VBSselection_Tight_Full},
        {"VBSBackgroundControl", VBSBackgroundControl},
        {"VBSBackgroundControlATLAS", VBSBackgroundControlATLAS},
        {"VBSBackgroundControl_Full", VBSBackgroundControl_Full},
        {"VBSBackgroundControlLoose", VBSBackgroundControlLoose},
        {"VBSBackgroundControlLoose_Full", VBSBackgroundControlLoose_Full},
        {"Inclusive2Jet", Inclusive2Jet},
        {"Inclusive2Jet_Full", Inclusive2Jet_Full},
        {"TightWithLooseVeto", TightWithLooseVeto}, 
        {"FourTopPlots", FourTopPlots},
        {"FourTopCutBasedEl", FourTopCutBasedEl},
	{"FourTopMVAEl", FourTopMVAEl}, 
    };

    std::map<std::string, Year> yearMap_ = {
        {"default", yrdefault},
        {"2016", yr2016},
        {"2017", yr2017},
        {"2018", yr2018},
    };
    
    std::map<std::string, Channel> channelMap_ = {
        {"e", e},                   {"m", m},         
        {"ep", ep},                 {"mp", mp},       {"ep", ep},       {"mn", mn},  
        {"ee", ee},                 {"em", em},       {"mm", mm},
        {"eee", eee},               {"eem", eem},     {"emm", emm},     {"mmm", mmm},
        {"eeee", eeee},             {"eemm", eemm},   {"mmee", mmee},   {"mmmm", mmmm},
        {"Inclusive", Inclusive},   {"lll", lll},     {"all", all},
    };


    std::vector<ChannelPair> allChannels_ = {};
    SystMap variations_ = {{Central, {}}};
    SystMap systematics_ = {};

    TList *currentHistDir_{nullptr};
    TH1D* sumWeightsHist_;

    std::vector<std::string> subprocesses_;
    bool doSystematics_;
    bool isNonprompt_ = false;
    bool applyScaleFactors_;
    bool applyPrefiringCorr_;
    
    // Readers to access the data (delete the ones you do not need).
    SelectorBase(TTree * /*tree*/ =0) { }
    virtual ~SelectorBase() { }
    virtual void    SetScaleFactors();
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    // We'll collect pointers to objects from derived classes
    // as they are registered with AddObject, and update them to
    // the new object when a dataset is switched
    std::set<TNamed**> allObjects_;
    // Derived classes override (and call) this to register new objects
    // With AddObject<Type>(localPtr, ...);
    virtual void SetupNewDirectory();
    void    SetBranches();
    void    LoadBranches(Long64_t entry, SystPair variation);
    virtual void    SetBranchesNanoAOD() {
        throw std::domain_error("NanoAOD ntuples not supported for selector!");
    }
    virtual void    LoadBranchesNanoAOD(Long64_t entry, SystPair variation) {
        throw std::domain_error("NanoAOD ntuples not supported for selector!");
    }
    virtual void    SetBranchesUWVV() {
        throw std::domain_error("UWVV ntuples not supported for selector!");
    }
    virtual void    LoadBranchesUWVV(Long64_t entry, SystPair variation) { 
        throw std::domain_error("UWVV ntuples not supported for selector!");
    }
    virtual void    SetBranchesBacon() {
        throw std::domain_error("Bacon ntuples not supported for selector!");
    }
    virtual void    LoadBranchesBacon(Long64_t entry, SystPair variation) {
        throw std::domain_error("Bacon ntuples not supported for selector!");
    }
    virtual void    FillHistograms(Long64_t entry, SystPair variation) { }
    void addSubprocesses(std::vector<std::string> processes);
    void makeOutputDirs();
    void setSubprocesses(std::string process);



    template<typename T, typename... Args>
	void AddObject(T* &ptr, Args... args) {
	static_assert(std::is_base_of<TNamed, T>::value, "Objects must inheirit from ROOT TNamed to be streamable from PROOF sessions");
	ptr = new T(args...);
	ptr->SetDirectory(gROOT);
	currentHistDir_->Add(ptr);
	allObjects_.insert((TNamed**) &ptr);
    };
    
    void UpdateDirectory();    
    ClassDef(SelectorBase,0);

 protected:
    // Maps to the histogram pointers themselves
    HistMap1D histMap1D_ = {};
    std::unordered_map<std::string, HistMap1D> subprocessHistMaps1D_ = {};
    std::vector<HistMap2D> subprocessWeightHistMaps1D_ = {};
    HistMap2D histMap2D_ = {};
    HistMap2D weighthistMap1D_ = {};
    HistMap3D weighthistMap2D_ {};

    std::vector<std::string> hists1D_ = {};
    std::vector<std::string> hists2D_ = {};
    std::vector<std::string> weighthists1D_ = {};
    std::vector<std::string> weighthists2D_ = {};
    // The histograms for which you also want systematic variations
    std::vector<std::string> systHists_ = {};
    std::vector<std::string> systHists2D_ = {};
    std::vector<Systematic> theoryVarSysts_ = {};

    // Variables
    std::string name_ = "Unnamed";
    std::string channelName_ = "Unnamed";
    Channel channel_ = Unknown;
    NtupleType ntupleType_ = NanoAOD;
    std::string selectionName_ = "tightleptons";
    Selection selection_ = tightleptons;
    Year year_ = yrdefault;
    bool isMC_;

    
    float GetPrefiringEfficiencyWeight(std::vector<float>* jetPt, std::vector<float>* jetEta);
    virtual std::string GetNameFromFile() { return ""; }
    void InitializeHistogramFromConfig(std::string name, ChannelPair channel, std::vector<std::string>& histData);
    void InitializeHistogramsFromConfig();
    std::vector<std::string> ReadHistDataFromConfig(std::string histDataString);
    std::string concatenateNames(const char* baseName, std::string& toAppend);
    std::string concatenateNames(const std::string& baseName, const char* toAppend);
    std::string concatenateNames(const std::string& baseName, std::string& toAppend);
    std::string getHistName(std::string histName, std::string variationName);
    std::string getHistName(std::string histName, std::string variationName, std::string channel);
    template<typename T>
    void InitializeHistMap(std::vector<std::string>& labels, std::unordered_map<HistLabel, T*>& histMap);

    // Filling Functions
    template<typename T, typename... Args>
	void SafeHistFill(std::unordered_map<HistLabel, T*>& container, 
		std::string name, Channel chan, Systematic var, Args... args) {
        SafeHistFill(container, name.c_str(), chan, var, args...);
    }

    template<typename T, typename... Args>
	void SafeHistFill(std::unordered_map<HistLabel, T*>& container, 
		const char* name, Channel chan, Systematic var, Args... args) {
        HistLabel histLabel = {name, chan, var};
        if (container[histLabel] != nullptr)
            container[histLabel]->Fill(args...);
    };
  
    template<typename T, typename... Args>
	void HistFullFill(std::unordered_map<HistLabel, T*>& container,
			  const char* histname, Systematic var, Args... args) {
	    SafeHistFill(container, histname, channel_, var, args...);
	    SafeHistFill(container, histname, all, var, args...);
    }
    template<typename T, typename... Args>
	void HistFullFill(std::unordered_map<HistLabel, T*>& container,
			  const std::string& histname, Systematic var, Args... args) {
        HistFullFill(container, histname.c_str(), var, args...);
    }
  
};

#endif

