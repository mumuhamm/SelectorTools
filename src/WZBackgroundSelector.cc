#include "Analysis/SelectorTools/interface/WZBackgroundSelector.h"
#include <TStyle.h>

void WZBackgroundSelector::SlaveBegin(TTree * /*tree*/)
{
    systematics_ = {
        {jetEnergyScaleUp, "CMS_scale_jUp"}, 
        {jetEnergyScaleDown, "CMS_scale_jDown"}, 
        {jetEnergyResolutionUp, "CMS_res_jUp"},
        {jetEnergyResolutionDown, "CMS_res_jDown"},
    };
    isNonpromptEstimate_ = true;
    doaQGC_ = false;
    WZSelector::SlaveBegin(0);
    fakeRate_allE_ = (ScaleFactor *) GetInputList()->FindObject("fakeRate_allE");
    if (fakeRate_allE_ == nullptr ) Abort("Must pass electron fake rate to input list!");
    fakeRate_allMu_ = (ScaleFactor *) GetInputList()->FindObject("fakeRate_allMu");
    if (fakeRate_allMu_ == nullptr ) Abort("Must pass muon fake rate to input list!");
}

void WZBackgroundSelector::SetupNewDirectory()
{
    WZSelector::SetupNewDirectory();
   
    // Insure that hist ranges are exactly the same as WZSelector, just change name
    for (const auto && obj : *currentHistDir_) {
        std::string name = obj->GetName();
        TNamed* named = dynamic_cast<TNamed*>(obj);
        named->SetName(name.insert(name.length()-3, "Fakes_").c_str());
    }
    
    AddObject<TH1D>(mjjHistPPF_, ("mjj_PPF_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistPFP_, ("mjj_PFP_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistFPP_, ("mjj_FPP_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistFFP_, ("mjj_FFP_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistFPF_, ("mjj_FPF_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistPFF_, ("mjj_PFF_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(mjjHistFFF_, ("mjj_FFF_"+channelName_).c_str(), "mjj; m_{jj} [GeV]; Events;", 15, 0, 1500);
    AddObject<TH1D>(ZMassHistPPF_, ("ZMass_PPF_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistPFP_, ("ZMass_PFP_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistFPP_, ("ZMass_FPP_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistFFP_, ("ZMass_FFP_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistFPF_, ("ZMass_FPF_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistPFF_, ("ZMass_PFF_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
    AddObject<TH1D>(ZMassHistFFF_, ("ZMass_FFF_"+channelName_).c_str(), "ZMass; m_{Z} [GeV]; Events;", 15, 76, 116);
}
float WZBackgroundSelector::getEventWeight() {
    float evtwgt = 0;
    if (!IsGenMatched3l())
        return true;

    if (IsFPPRegion()) {
        if (!isVBS_ || PassesVBSSelection(true))
            mjjHistFPP_->Fill(mjj, weight);
        evtwgt = getl1FakeRate()*weight;
    }
    else if (IsPFPRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistPFP_->Fill(mjj, weight);
            ZMassHistPFP_->Fill(ZMass, weight);
        }
        evtwgt = getl2FakeRate()*weight;
    }
    else if (IsPPFRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistPPF_->Fill(mjj, weight);
            ZMassHistPPF_->Fill(ZMass, weight);
        }
        evtwgt = getl3FakeRate()*weight;
    }
    else if (IsFFFRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistFFF_->Fill(mjj, weight);
            ZMassHistFFF_->Fill(ZMass, weight);
        }
        evtwgt = getl1FakeRate()*getl2FakeRate()*getl3FakeRate()*weight;
    }
    else if (IsFPFRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistFPF_->Fill(mjj, weight);
            ZMassHistFPF_->Fill(ZMass, weight);
        }
        evtwgt = -1*getl1FakeRate()*getl3FakeRate()*weight;
    }
    else if (IsFFPRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistFFP_->Fill(mjj, weight);
            ZMassHistFFP_->Fill(ZMass, weight);
        }
        evtwgt = -1*getl1FakeRate()*getl2FakeRate()*weight;
    }
    else if (IsPFFRegion()) {
        if (!isVBS_ || PassesVBSSelection(true)) {
            mjjHistPFF_->Fill(mjj, weight);
            ZMassHistPFF_->Fill(ZMass, weight);
        }
        evtwgt = -1*getl2FakeRate()*getl3FakeRate()*weight;
    }
    if (isMC_)
        evtwgt *= GetPrefiringEfficiencyWeight(jetPt, jetEta);
    return evtwgt;
}

void WZBackgroundSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {
    WZSelector::LoadBranches(entry, variation);
    weight = getEventWeight();
}

float WZBackgroundSelector::getl1FakeRate() {
    float pt_fillval = l1Pt < FR_MAX_PT_ ? l1Pt : FR_MAX_PT_ - 0.01;
    float fr = 1;
    if (channel_ == eem|| channel_ == eee)
        fr = fakeRate_allE_->Evaluate2D(pt_fillval, std::abs(l1Eta));
    else
        fr = fakeRate_allMu_->Evaluate2D(pt_fillval, std::abs(l1Eta));
    return fr/(1-fr); 
}

float WZBackgroundSelector::getl2FakeRate() {
    float pt_fillval = l2Pt < FR_MAX_PT_ ? l2Pt : FR_MAX_PT_ - 0.01;
    float fr = 1;
    if (channel_ == eem || channel_ == eee)
        fr = fakeRate_allE_->Evaluate2D(pt_fillval, std::abs(l2Eta));
    else
        fr = fakeRate_allMu_->Evaluate2D(pt_fillval, std::abs(l2Eta));
    return fr/(1-fr); 
}

float WZBackgroundSelector::getl3FakeRate() {
    float pt_fillval = l3Pt < FR_MAX_PT_ ? l3Pt : FR_MAX_PT_ - 0.01;
    float fr = 1;
    if (channel_ == emm || channel_ == eee)
        fr = fakeRate_allE_->Evaluate2D(pt_fillval, std::abs(l3Eta));
    else
        fr = fakeRate_allMu_->Evaluate2D(pt_fillval, std::abs(l3Eta));
    return fr/(1-fr); 
}

bool WZBackgroundSelector::IsFPPRegion() {
    return !zlep1IsTight() && zlep2IsTight() && lepton3IsTight();
}

bool WZBackgroundSelector::IsPFPRegion() {
    return zlep1IsTight() && !zlep2IsTight() && lepton3IsTight();
}

bool WZBackgroundSelector::IsPPFRegion() {
    return zlep1IsTight() && zlep2IsTight() && !lepton3IsTight();
}

bool WZBackgroundSelector::IsFPFRegion() {
    return !zlep1IsTight() && zlep2IsTight() && !lepton3IsTight();
}

bool WZBackgroundSelector::IsFFPRegion() {
    return !zlep1IsTight() && !zlep2IsTight() && lepton3IsTight();
}

bool WZBackgroundSelector::IsPFFRegion() {
    return zlep1IsTight() && !zlep2IsTight() && !lepton3IsTight();
}

bool WZBackgroundSelector::IsFFFRegion() {
    return !zlep1IsTight() && !zlep2IsTight() && !lepton3IsTight();
}
