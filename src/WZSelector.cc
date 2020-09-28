#include "Analysis/SelectorTools/interface/WZSelector.h"
#include "TLorentzVector.h"

void WZSelector::Init(TTree *tree)
{
    weight_info_ = 0;
    if (isMC_) {
        weight_info_ = GetLheWeightInfo();
    }

    systematics_ = {
        {jetEnergyScaleUp, "CMS_scale_jUp"}, 
        {jetEnergyScaleDown, "CMS_scale_jDown"}, 
        {jetEnergyResolutionUp, "CMS_res_jUp"},
        {jetEnergyResolutionDown, "CMS_res_jDown"},
        {metUnclusteredEnergyUp, "CMS_scale_unclEnergyUp"},
        {metUnclusteredEnergyDown, "CMS_scale_unclEnergyDown"},
        {muonScaleUp, "CMS_scale_mUp"},
        {muonScaleDown, "CMS_scale_mDown"},
        {electronScaleUp, "CMS_scale_eUp"},
        {electronScaleDown, "CMS_scale_eDown"},
        {pileupUp, "CMS_pileupUp"},
        {pileupDown, "CMS_pileupDown"},
        {electronEfficiencyUp, "CMS_eff_eUp"},
        {electronEfficiencyDown, "CMS_eff_eDown"},
        {muonEfficiencyUp, "CMS_eff_mUp"},
        {muonEfficiencyDown, "CMS_eff_mDown"},
    };

    systHists_ = {
        "yield",
        //"backgroundControlYield",
        //"mjj",
        //"dEtajj",
        //"MTWZ",
    };

    hists1D_ = {"yield", "Zlep1_Eta", "Zlep1_Phi", "Zlep1_Pt",
        "Zlep2_Eta", "Zlep2_Phi", "Zlep2_Pt", "Wlep_Eta", "Wlep_Phi", "Wlep_Pt",
        "ZMass", "Mass", "MET", "nJets",
        //"ZPt",
        //"ZEta",
        //"ZPhi",
        //"dR_lW_Z",
        //"backgroundControlYield",
        //"nTruePU",
        //"nvtx",
        //"MTWZ",
        //"M3lMET",
        //"Mass",
        //"Pt",
        //"nJets",
        //"nJetCSVv2T",
        //"jetPt[0]",
        //"jetPt[1]",
        //"jetPt[2]",
        //"jetEta12",
        //"jetEta[0]",
        //"jetEta[1]",
        //"jetEta[2]",
        //"mjj",
        //"MtW",
        //"dEtajj",
        //"dRjj",
        //"zep3l",
        //"zepj3",
        //"Eta",
        //"m_l1l3",
    };

    WZSelectorBase::Init(tree);
}

// Values from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceScaleResolRun2
float WZSelector::GetMuonScaleUncertainty(float muEta) {
    if (muEta < -2.1)
        return 0.027;
    if (std::abs(muEta) < 2.1 && std::abs(muEta) > 1.2)
        return 0.009;
    if (std::abs(muEta) < 1.2)
        return 0.004;
    if (muEta > 2.1)
        return 0.0017;
    throw std::out_of_range("Muon eta out of range of possible values! Range was " + std::to_string(muEta));
}

void WZSelector::SetBranchesUWVV() {
    WZSelectorBase::SetBranchesUWVV();
    if (isMC_) {
        if (weight_info_ > 0)
            fChain->SetBranchAddress("scaleWeights", &scaleWeights, &b_scaleWeights);
        if ((weight_info_ == 2 || weight_info_ == 3) && doSystematics_)
            fChain->SetBranchAddress("pdfWeights", &pdfWeights, &b_pdfWeights);
        fChain->SetBranchAddress("mjj_jesUp", &mjj_jesUp, &b_mjj_jesUp);
        fChain->SetBranchAddress("mjj_jesDown", &mjj_jesDown, &b_mjj_jesDown);
        fChain->SetBranchAddress("mjj_jerUp", &mjj_jerUp, &b_mjj_jerUp);
        fChain->SetBranchAddress("mjj_jerDown", &mjj_jerDown, &b_mjj_jerDown);
        fChain->SetBranchAddress("jetEta_jesUp", &jetEta_jesUp, &b_jetEta_jesUp);
        fChain->SetBranchAddress("jetEta_jesDown", &jetEta_jesDown, &b_jetEta_jesDown);
        fChain->SetBranchAddress("jetEta_jerUp", &jetEta_jerUp, &b_jetEta_jerUp);
        fChain->SetBranchAddress("jetEta_jerDown", &jetEta_jerDown, &b_jetEta_jerDown);
        fChain->SetBranchAddress("jetPt_jesUp", &jetPt_jesUp, &b_jetPt_jesUp);
        fChain->SetBranchAddress("jetPt_jesDown", &jetPt_jesDown, &b_jetPt_jesDown);
        fChain->SetBranchAddress("jetPt_jerUp", &jetPt_jerUp, &b_jetPt_jerUp);
        fChain->SetBranchAddress("jetPt_jerDown", &jetPt_jerDown, &b_jetPt_jerDown);
        fChain->SetBranchAddress("type1_pfMETEt_jesUp", &type1_pfMETEt_jesUp, &b_type1_pfMETEt_jesUp);
        fChain->SetBranchAddress("type1_pfMETEt_jesDown", &type1_pfMETEt_jesDown, &b_type1_pfMETEt_jesDown);
        fChain->SetBranchAddress("type1_pfMETEt_jerUp", &type1_pfMETEt_jerUp, &b_type1_pfMETEt_jerUp);
        fChain->SetBranchAddress("type1_pfMETEt_jerDown", &type1_pfMETEt_jerDown, &b_type1_pfMETEt_jerDown);
        fChain->SetBranchAddress("type1_pfMETEt_unclusteredEnUp", &type1_pfMETEt_unclusteredEnUp, &b_type1_pfMETEt_unclusteredEnUp);
        fChain->SetBranchAddress("type1_pfMETEt_unclusteredEnDown", &type1_pfMETEt_unclusteredEnDown, &b_type1_pfMETEt_unclusteredEnDown);
        fChain->SetBranchAddress("type1_pfMETEt_UncTool", &type1_pfMETEt_UncTool, &b_type1_pfMETEt_UncTool);
        if (channel_ == eee) {
            fChain->SetBranchAddress("e1scaleCorrError", &l1PtScaleCorrErr, &b_l1PtScaleCorrErr);
            fChain->SetBranchAddress("e2scaleCorrError", &l2PtScaleCorrErr, &b_l2PtScaleCorrErr);
            fChain->SetBranchAddress("e3scaleCorrError", &l3PtScaleCorrErr, &b_l3PtScaleCorrErr);
        }
        else if (channel_ == eem) {
            fChain->SetBranchAddress("e1scaleCorrError", &l1PtScaleCorrErr, &b_l1PtScaleCorrErr);
            fChain->SetBranchAddress("e2scaleCorrError", &l2PtScaleCorrErr, &b_l2PtScaleCorrErr);
        }
        else if (channel_ == emm) {
            fChain->SetBranchAddress("escaleCorrError", &l3PtScaleCorrErr, &b_l3PtScaleCorrErr);
        }
    }
    
    fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
    fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
    fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
    fChain->SetBranchAddress("jetCSVv2", &jetCSVv2, &b_jetCSVv2);
    fChain->SetBranchAddress("Eta", &Eta, &b_Eta);
    fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
    fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
    fChain->SetBranchAddress("mjj", &mjj, &b_mjj);

    if (channel_ == eee) {
        fChain->SetBranchAddress("e1_e3_Mass", &Zlep1_Wlep_Mass, &b_Zlep1_Wlep_Mass);
        fChain->SetBranchAddress("e2_e3_Mass", &Zlep2_Wlep_Mass, &b_Zlep2_Wlep_Mass);
        fChain->SetBranchAddress("e1_e2_Pt", &ZPt, &b_ZPt);
        fChain->SetBranchAddress("e1_e2_Eta", &ZEta, &b_ZEta);
        fChain->SetBranchAddress("e1_e2_Phi", &ZPhi, &b_ZPhi);
    }
    else if (channel_ == eem) { 
        fChain->SetBranchAddress("e1_m_Mass", &Zlep1_Wlep_Mass, &b_Zlep1_Wlep_Mass);
        fChain->SetBranchAddress("e2_m_Mass", &Zlep2_Wlep_Mass, &b_Zlep2_Wlep_Mass);
        fChain->SetBranchAddress("e1_e2_Pt", &ZPt, &b_ZPt);
        fChain->SetBranchAddress("e1_e2_Eta", &ZEta, &b_ZEta);
        fChain->SetBranchAddress("e1_e2_Phi", &ZPhi, &b_ZPhi);
    }
    else if (channel_ == emm) { 
        fChain->SetBranchAddress("e_m1_Mass", &Zlep1_Wlep_Mass, &b_Zlep1_Wlep_Mass);
        fChain->SetBranchAddress("e_m2_Mass", &Zlep2_Wlep_Mass, &b_Zlep2_Wlep_Mass);
        fChain->SetBranchAddress("m1_m2_Pt", &ZPt, &b_ZPt);
        fChain->SetBranchAddress("m1_m2_Eta", &ZEta, &b_ZEta);
        fChain->SetBranchAddress("m1_m2_Phi", &ZPhi, &b_ZPhi);
    }
    else if (channel_ == mmm) { 
        fChain->SetBranchAddress("m1_m3_Mass", &Zlep1_Wlep_Mass, &b_Zlep1_Wlep_Mass);
        fChain->SetBranchAddress("m2_m3_Mass", &Zlep2_Wlep_Mass, &b_Zlep2_Wlep_Mass);
        fChain->SetBranchAddress("m1_m2_Pt", &ZPt, &b_ZPt);
        fChain->SetBranchAddress("m1_m2_Eta", &ZEta, &b_ZEta);
        fChain->SetBranchAddress("m1_m2_Phi", &ZPhi, &b_ZPhi);
    }

}

unsigned int WZSelector::GetLheWeightInfo() {
    std::vector<std::string> noLheWeights = {
        "ggZZ2e2mu", "ggZZ4e", "ggZZ4m", "wzjj-vbfnlo-of",
        "wzjj-vbfnlo-sf", "st-schan", "st-tchan", "st-tchan-t", "st-tchan-tbar"
    };
    std::vector<std::string> scaleAndPdfWeights = {
        "wz3lnu-powheg", "wz3lnu-mg5amcnlo",
        "wz3lnu-mgmlm-0j", "wz3lnu-mgmlm-1j",
        "wz3lnu-mgmlm-2j", "wz3lnu-mgmlm-3j", "wlljj-ewk", 
        "tzq", "ww", "www", "wwz", "zz-powheg",
        "zz4l-mg5amcnlo", "zz4ljj-ewk", "zz2l2vjj-ewk"
    };
    std::vector<std::string> allLheWeights = {
        "wzjj-aqgcft", "wzjj-aqgcfm", "wzjj-aqgcfs",
        //"wz-atgc_pt0-200", "wz-atgc_pt200-300",
        //"wz-atgc_pt300"
    };

    if (std::find(noLheWeights.begin(), noLheWeights.end(), name_) != noLheWeights.end())
        return 0;
    if (std::find(scaleAndPdfWeights.begin(), scaleAndPdfWeights.end(), name_) != scaleAndPdfWeights.end())
        return 2;
    if (std::find(allLheWeights.begin(), allLheWeights.end(), name_) != allLheWeights.end())
        return 3;
    return 1;
}

void WZSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    WZSelectorBase::LoadBranchesUWVV(entry, variation);
    //weight *= GetPrefiringEfficiencyWeight(jetPt, jetEta);

    //b_MtToMET->GetEntry(entry);
    b_ZPhi->GetEntry(entry);
    b_ZEta->GetEntry(entry);
    b_jetPt->GetEntry(entry);
    b_jetEta->GetEntry(entry);
    b_jetPhi->GetEntry(entry);
    b_Eta->GetEntry(entry);
    b_mjj->GetEntry(entry);
    if (variation.first == Central) {
        if (isMC_ && doSystematics_) {
            if (isMC_ && weight_info_ > 0) {
                b_scaleWeights->GetEntry(entry);
                lheWeights = *scaleWeights;
                if (weight_info_ == 2) {
                    b_pdfWeights->GetEntry(entry);
                    // Only keep NNPDF weights
                    lheWeights.insert(lheWeights.end(), pdfWeights->begin(), 
                        pdfWeights->begin()+std::min(static_cast<size_t>(103), pdfWeights->size()));
                }
                else if (weight_info_ == 3) {
                    b_pdfWeights->GetEntry(entry);
                    lheWeights.insert(lheWeights.end(), pdfWeights->begin(), pdfWeights->end());
                }
            }
        }

        if (histMap1D_[{"MTWZ", channel_, Central}] != nullptr || histMap1D_[{"M3lMET", channel_, Central}] == nullptr) {
            TLorentzVector l1 = TLorentzVector();
            l1.SetPtEtaPhiM(l1Pt, l1Eta, l1Phi, 0);
            TLorentzVector l2 = TLorentzVector();
            l2.SetPtEtaPhiM(l2Pt, l2Eta, l2Phi, 0);
            TLorentzVector l3 = TLorentzVector();
            l3.SetPtEtaPhiM(l3Pt, l3Eta, l3Phi, 0);
            TLorentzVector met = TLorentzVector();
            met.SetPtEtaPhiM(MET, 0, type1_pfMETPhi, 0);
            TLorentzVector zp4 = l1+l2;
            TLorentzVector wp4 = l3+met;
            TLorentzVector wzp4 = l1+l2+l3+met;
            float sumEt = zp4.Et() + wp4.Et();
            MtWZ = std::sqrt(sumEt*sumEt - wzp4.Pt()*wzp4.Pt());
            M3lMET = (l1+l2+l3+met).M();
        }
    }
    else if (isMC_) {
        if (variation.first == jetEnergyScaleUp) {
            b_mjj_jesUp->GetEntry(entry);
            b_jetEta_jesUp->GetEntry(entry);
            b_jetPt_jesUp->GetEntry(entry);
            b_type1_pfMETEt_jesUp->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);

            mjj = mjj_jesUp;
            jetEta = jetEta_jesUp;
            jetPt = jetPt_jesUp;
            MET = type1_pfMETEt_jesUp*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == jetEnergyScaleDown) {
            b_mjj_jesDown->GetEntry(entry);
            b_jetEta_jesDown->GetEntry(entry);
            b_jetPt_jesDown->GetEntry(entry);
            b_type1_pfMETEt_jesDown->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);

            mjj = mjj_jesDown;
            jetEta = jetEta_jesDown;
            jetPt = jetPt_jesDown;
            MET = type1_pfMETEt_jesDown*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == jetEnergyResolutionUp) {
            b_mjj_jerUp->GetEntry(entry);
            b_jetEta_jerUp->GetEntry(entry);
            b_jetPt_jerUp->GetEntry(entry);
            b_type1_pfMETEt_jerUp->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);

            mjj = mjj_jerUp;
            jetEta = jetEta_jerUp;
            jetPt = jetPt_jerUp;
            MET = type1_pfMETEt_jerUp*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == jetEnergyResolutionDown) {
            b_mjj_jerDown->GetEntry(entry);
            b_jetEta_jerDown->GetEntry(entry);
            b_jetPt_jerDown->GetEntry(entry);
            b_type1_pfMETEt_jerDown->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);

            mjj = mjj_jerDown;
            jetEta = jetEta_jerDown;
            jetPt = jetPt_jerDown;
            MET = type1_pfMETEt_jerDown*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == muonScaleUp || variation.first == muonScaleDown) {
            if (channel_ == eee)
                return;
            if (channel_ == eem) {
                float scaleUnc = GetMuonScaleUncertainty(l3Eta);
                l3Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
            }
            else if (channel_ == emm) {
                float scaleUnc = GetMuonScaleUncertainty(l1Eta);
                l1Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
                scaleUnc = GetMuonScaleUncertainty(l2Eta);
                l2Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
            }
            else if (channel_ == mmm) {
                float scaleUnc = GetMuonScaleUncertainty(l1Eta);
                l1Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
                scaleUnc = GetMuonScaleUncertainty(l2Eta);
                l2Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
                scaleUnc = GetMuonScaleUncertainty(l3Eta);
                l3Pt *= variation.first == muonScaleUp ? (1+scaleUnc) : (1-scaleUnc);
            }
            SetMasses();
        }
        else if (variation.first == electronScaleUp || variation.first == electronScaleDown) {
            if (channel_ == eee) {
                b_l1PtScaleCorrErr->GetEntry(entry);
                b_l2PtScaleCorrErr->GetEntry(entry);
                b_l3PtScaleCorrErr->GetEntry(entry);
                l1Pt *= variation.first == electronScaleUp ? (1+l1PtScaleCorrErr) : (1-l1PtScaleCorrErr);
                l2Pt *= variation.first == electronScaleUp ? (1+l2PtScaleCorrErr) : (1-l2PtScaleCorrErr);
                l3Pt *= variation.first == electronScaleUp ? (1+l3PtScaleCorrErr) : (1-l3PtScaleCorrErr);
            }
            if (channel_ == eem) {
                b_l1PtScaleCorrErr->GetEntry(entry);
                b_l2PtScaleCorrErr->GetEntry(entry);
                l1Pt *= variation.first == electronScaleUp ? (1+l1PtScaleCorrErr) : (1-l1PtScaleCorrErr);
                l2Pt *= variation.first == electronScaleUp ? (1+l2PtScaleCorrErr) : (1-l2PtScaleCorrErr);
            }
            else if (channel_ == emm) {
                b_l3PtScaleCorrErr->GetEntry(entry);
                l3Pt *= variation.first == electronScaleUp ? (1+l3PtScaleCorrErr) : (1-l3PtScaleCorrErr);
            }
            else if (channel_ == mmm)
                return;
            SetMasses();
        }
        else if (variation.first == metUnclusteredEnergyDown) {
            b_type1_pfMETEt_unclusteredEnUp->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);
            MET = type1_pfMETEt_unclusteredEnUp*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == metUnclusteredEnergyDown) {
            b_type1_pfMETEt_unclusteredEnDown->GetEntry(entry);
            b_type1_pfMETEt_UncTool->GetEntry(entry);
            MET = type1_pfMETEt_unclusteredEnDown*MET/type1_pfMETEt_UncTool;
        }
        else if (variation.first == pileupUp) {
            weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftUp)/pileupSF_->Evaluate1D(nTruePU);
        }
        else if (variation.first == pileupDown) {
            weight *= pileupSF_->Evaluate1D(nTruePU, ScaleFactor::ShiftDown)/pileupSF_->Evaluate1D(nTruePU);
        }
        else if (variation.first == electronEfficiencyUp || variation.first == electronEfficiencyDown ||
                    variation.first == muonEfficiencyUp || variation.first == muonEfficiencyDown) {
            ShiftEfficiencies(variation.first);
        }
    }

    auto deltaRjj = [](std::vector<float>* jEta, std::vector<float>* jPhi) {
        if (jEta->size() != jPhi->size() || jEta->size() < 2)
            return -1.;
        double etaDiff = jEta->at(0) - jEta->at(1);
        double phiDiff = jPhi->at(0) - jPhi->at(1);
        return std::sqrt(etaDiff*etaDiff + phiDiff*phiDiff);
    };

    auto deltaEtajj = [](std::vector<float>* jEta) {
        if (jEta->size() < 2)
            return -1.;
        double etaDiff = jEta->at(0) - jEta->at(1);
        return std::abs(etaDiff);
    };

    auto zeppenfeld = [](std::vector<float>* jEta, float objEta) {
        if (jEta->size() < 2)
            return -999.;
        return objEta - 0.5*(jEta->at(1) + jEta->at(0));
    };


    dEtajj = deltaEtajj(jetEta);
    dRjj = deltaRjj(jetEta, jetPhi);
    zep3l = zeppenfeld(jetEta, Eta);
}

bool WZSelector::PassesVBSBackgroundControlSelection() {
    if (jetPt->size() != jetEta->size() || jetPt->size() < 2)
        return false;

    if (selection_ != VBSselection_Loose &&
            selection_ != VBSselection_Loose_Full) {
        if (jetPt->at(0) < 50 || jetPt->at(1) < 50)
            return false;
    }
    else
        return (mjj > 100 && (mjj < 500 || dEtajj < 2.5));

    return (mjj > 100 && (mjj < 500 || dEtajj < 2.5 || std::abs(zep3l) > 2.5));
}

void WZSelector::ShiftEfficiencies(Systematic variation) {
    ScaleFactor::Variation shift = ScaleFactor::Variation::ShiftUp;
    if (variation == electronEfficiencyDown || variation == muonEfficiencyDown)
        shift = ScaleFactor::Variation::ShiftDown;

    if (channel_ == eee && (variation == electronEfficiencyUp || variation == electronEfficiencyDown)) {
        weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/eIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        weight *= eGsfSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    }
    else if (channel_ == eem) {
        if (variation == electronEfficiencyUp || variation == electronEfficiencyDown) {
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/eIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        }
        else if (variation == muonEfficiencyUp || variation == muonEfficiencyDown) {
            weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
    }
    else if (channel_ == emm) {
        if (variation == electronEfficiencyUp || variation == electronEfficiencyDown) {
            weight *= eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
            weight *= eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/eGsfSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        }
        else if (variation == muonEfficiencyUp || variation == muonEfficiencyDown) {
            weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
            weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
            weight *= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        }
    }
    else if (channel_ == mmm && (variation == muonEfficiencyUp || variation == muonEfficiencyDown)) {
        weight *= mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/mIdSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l1Eta), l1Pt);
        weight *= mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/mIdSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l2Eta), l2Pt);
        weight *= mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIdSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
        weight *= mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt, shift)/mIsoSF_->Evaluate2D(std::abs(l3Eta), l3Pt);
    }
}

bool WZSelector::PassesVBSSelection(bool noBlind) {
    if (!(isMC_ || noBlind))
        return false;
    if (jetPt->size() != jetEta->size() || jetPt->size() < 2)
        return false;

    // Use optimized point of pT(j1,j2) > 50 GeV
    if (selection_ != VBSselection_Loose && 
            selection_ != VBSselection_Loose_Full &&
            selection_ != VBSBackgroundControl) { // &&
        if (jetPt->at(0) < 50 || jetPt->at(1) < 50)
            return false;
    }
    else if (jetPt->at(0) < 40 || jetPt->at(1) < 40)
        return false;
    if (selection_ == VBSselection_Tight ||
                selection_ == VBSselection_Tight_Full) { 
        if (std::abs(zep3l) > 2.5)
            return false;
        return mjj > 500 && dEtajj > 2.5;
    }
    // Background control
    else if (selection_ == VBSBackgroundControl) { 
        return ((mjj > 500 && dEtajj < 2.5) ||
                (mjj < 500 && dEtajj > 2.5));
    }
    else if (selection_ == VBSBackgroundControlATLAS) { 
        return jetPt->at(0) > 40 && jetPt->at(1) > 40 && mjj > 150 && mjj < 500;
    }
    else if (selection_ == VBSBackgroundControlLoose ||
             selection_ == VBSBackgroundControlLoose_Full) { 
        return PassesVBSBackgroundControlSelection();
    }
    return mjj > 500 && dEtajj > 2.5;
    // ATLAS
    //return mjj > 150;
}

bool WZSelector::PassesFullWZSelection(Long64_t entry) {
    if (ZMass > 106.1876 || ZMass < 76.1876)
        return false;
    if (l1Pt < 25 || l2Pt < 15 || l3Pt < 20)
        return false;
    if (Mass < 100)
        return false;
    if (MET < 30)
        return false;

    //b_jetCSVv2->GetEntry(entry);
    //for (const auto& jetCSVval : *jetCSVv2) {
    //    if (jetCSVval > 0.9535)
    //        return false;
    //}
    //b_Zlep1_Wlep_Mass->GetEntry(entry);
    //b_Zlep2_Wlep_Mass->GetEntry(entry);
    //if (Zlep1_Wlep_Mass < 4 || Zlep2_Wlep_Mass < 4)
    //    return false;

    return true;
}

bool WZSelector::PassesBaseSelection(Long64_t entry, bool tightLeps, Selection selection) { 
    //if (!(Flag_BadChargedCandidateFilterPass
    //        && Flag_HBHENoiseFilterPass 
    //        && Flag_HBHENoiseIsoFilterPass 
    //        && Flag_BadPFMuonFilterPass
    //        && Flag_EcalDeadCellTriggerPrimitiveFilterPass 
    //        && Flag_goodVerticesPass 
    //        && Flag_globalTightHalo2016FilterPass
    //        && (isMC_ || Flag_eeBadScFilterPass
    //            //&& !Flag_duplicateMuonsPass 
    //            //&& !Flag_badMuonsPass)
    //            // No longer vetoing events failing these filters
    //            // (we trust the MET in ReMiniAOD)
    //        )
    //    )
    //) 
    //    return false;
    if (!passesLeptonVeto)
        return false;
    if (applyFullSelection_ && !PassesFullWZSelection(entry))
        return false;
    else if (selection == FakeRateSelectionLoose || selection == FakeRateSelectionTight) {
        if (l1Pt < 25 || l2Pt < 15)
            return false;
        if (ZMass > 111.1876 || ZMass < 81.1876)
            return false;
        if (MET > 25)
            return false;
        if (l3MtToMET > 30)
            return false;
        if (!tightZLeptons())
            return false;
        if (selection == FakeRateSelectionLoose)
            tightLeps = false;
    }
    else if ((selection == Inclusive2Jet || selection == Inclusive2Jet_Full) &&
                (jetPt->size() < 2 || jetPt->at(0) < 50 || jetPt->at(1) < 50))
        return false;
    
    if (tightLeps && !(zlep1IsTight() && zlep2IsTight() && lepton3IsTight()))
        return false;
    //if (!IsGenMatched3l()) {
    //    return false;
    //}

    return true;
}

void WZSelector::FillVBSHistograms(float weight, bool noBlind,
        std::pair<Systematic, std::string> variation) { 
    // JES/JER uncertainties
    // Need to separate check VBS cuts using JER/JES variations
    SafeHistFill(histMap2D_, "mjj_etajj_2D", channel_, variation.first, 
        mjj, dEtajj, weight*(isMC_ || noBlind || mjj < 500 || dEtajj < 2.5));
    SafeHistFill(histMap1D_, "zep3l", channel_, variation.first, zep3l, weight);
    SafeHistFill(histMap2D_, "mjj_dRjj_2D", channel_, variation.first, mjj, dRjj, weight*(isMC_ || noBlind || mjj < 500 || dEtajj < 2.5));

    SafeHistFill(histMap1D_, "mjj", channel_, variation.first, mjj, weight*(isMC_ || (mjj < 500) || noBlind));
    SafeHistFill(histMap1D_, "dEtajj", channel_, variation.first, dEtajj, weight*(isMC_ || (dEtajj < 2.5) || noBlind));
    SafeHistFill(histMap1D_, "dRjj", channel_, variation.first, dRjj, weight*(isMC_ || (dRjj < 2.5) || noBlind));

    if (jetPt->size() > 0 && jetPt->size() == jetEta->size()) {
        SafeHistFill(histMap1D_, "jetPt[0]", channel_, variation.first, jetPt->at(0), weight);
        SafeHistFill(histMap1D_, "jetEta[0]", channel_, variation.first, jetEta->at(0), weight);
    }
    if (jetPt->size() > 1 && jetPt->size() == jetEta->size()) {
        SafeHistFill(histMap1D_, "jetPt[1]", channel_, variation.first, jetPt->at(1), weight);
        SafeHistFill(histMap1D_, "jetEta[1]", channel_, variation.first, jetEta->at(1), weight);
    }
    if (jetPt->size() > 2 && jetPt->size() == jetEta->size()) {
        SafeHistFill(histMap1D_, "jetPt[2]", channel_, variation.first, jetPt->at(2), weight);
        SafeHistFill(histMap1D_, "jetEta[2]", channel_, variation.first, jetEta->at(2), weight);
    }
     
    if (jetEta->size() > 3)
        SafeHistFill(histMap1D_, "zepj3", channel_, variation.first, jetEta->at(2) - 0.5*(jetEta->at(1) + jetEta->at(0)), weight);
    
    if (histMap1D_[{"jetEta12", channel_, variation.first}] != nullptr && jetEta->size() > 1) {
        histMap1D_[{"jetEta12", channel_, variation.first}]->Fill(jetEta->at(0), weight);
        histMap1D_[{"jetEta12", channel_, variation.first}]->Fill(jetEta->at(1), weight);
    }
}

void WZSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    bool noBlind = true;
    if (!PassesBaseSelection(entry, true, selection_))
        return;
    //bool passesVBS = PassesVBSSelection(noBlind);
    bool passesVBS = false;

    //if (histMap1D_["backgroundControlYield", channel_, variation.first] != nullptr)
    //    if (PassesVBSBackgroundControlSelection())
    //        histMap1D_["backgroundControlYield", channel_, variation.first]->Fill(1, weight);

    if ((variation.first == Central || (doaQGC_ && isaQGC_)) && isMC_) {
        for (size_t i = 0; i < lheWeights.size(); i++) {
            if (PassesVBSBackgroundControlSelection())
                SafeHistFill(weighthistMap1D_, "backgroundControlYield", channel_, variation.first, 1, i, lheWeights[i]/lheWeights[0]*weight);
            if (isVBS_ && !passesVBS)
                continue;
            SafeHistFill(weighthistMap1D_, "yield", channel_, variation.first, 1, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "mjj", channel_, variation.first, mjj, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "m_l1l3", channel_, variation.first, Zlep1_Wlep_Mass, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "m_l2l3", channel_, variation.first, Zlep2_Wlep_Mass, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "Mass", channel_, variation.first, Mass, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "MTWZ", channel_, variation.first, MtWZ, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "M3lMET", channel_, variation.first, M3lMET, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "l1Pt", channel_, variation.first, l1Pt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "l2Pt", channel_, variation.first, l2Pt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "l3Pt", channel_, variation.first, l3Pt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "ZPt", channel_, variation.first, ZPt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "Pt", channel_, variation.first, ZPt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap1D_, "Mass", channel_, variation.first, ZPt, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap2D_, "mjj_etajj_2D", channel_, variation.first, mjj, dEtajj, i, lheWeights[i]/lheWeights[0]*weight);
            SafeHistFill(weighthistMap2D_, "mjj_dRjj_2D", channel_, variation.first, mjj, dRjj, i, lheWeights[i]/lheWeights[0]*weight);
        }
    }
    if (isVBS_ && !passesVBS)
        return;
    //FillVBSHistograms(weight, noBlind, variation);

    SafeHistFill(histMap1D_, "yield", channel_, variation.first, 1, weight);
    SafeHistFill(histMap1D_, "Mass", channel_, variation.first, Mass, 
        weight*(isMC_ || Mass < 400 || noBlind));
    SafeHistFill(histMap1D_, "ZMass", channel_, variation.first, ZMass, weight);
    SafeHistFill(histMap1D_, "Zlep1_Pt", channel_, variation.first, l1Pt, weight);
    SafeHistFill(histMap1D_, "Zlep1_Eta", channel_, variation.first, l1Eta, weight);
    SafeHistFill(histMap1D_, "Zlep1_Phi", channel_, variation.first, l1Phi, weight);
    SafeHistFill(histMap1D_, "Zlep2_Pt", channel_, variation.first, l2Pt, weight);
    SafeHistFill(histMap1D_, "Zlep2_Eta", channel_, variation.first, l2Eta, weight);
    SafeHistFill(histMap1D_, "Zlep2_Phi", channel_, variation.first, l2Phi, weight);
    SafeHistFill(histMap1D_, "Wlep_Pt", channel_, variation.first, l3Pt, weight);
    SafeHistFill(histMap1D_, "Wlep_Eta", channel_, variation.first, l3Eta, weight);
    SafeHistFill(histMap1D_, "Wlep_Phi", channel_, variation.first, l3Phi, weight);
    SafeHistFill(histMap1D_, "MET", channel_, variation.first, MET, weight);
    // Just doing what works for now
    return;
    SafeHistFill(histMap1D_, "nJets", channel_, variation.first, jetPt->size(), weight);
    SafeHistFill(histMap1D_, "m_l1l3", channel_, variation.first, Zlep1_Wlep_Mass, weight);
    SafeHistFill(histMap1D_, "m_l2l3", channel_, variation.first, Zlep2_Wlep_Mass, weight);
    SafeHistFill(histMap1D_, "ZPhi", channel_, variation.first, ZPhi, weight);
    SafeHistFill(histMap1D_, "ZEta", channel_, variation.first, ZEta, weight);
    SafeHistFill(histMap1D_, "MtW", channel_, variation.first, l3MtToMET, weight);
    SafeHistFill(histMap1D_, "Eta", channel_, variation.first, Eta, weight);

    if (histMap1D_[{"dR_lW_Z", channel_, variation.first}] != nullptr) {
        float dPhi_lW_Z = ZPhi - l3Phi;
        float dEta_lW_Z = ZEta - l3Eta;
        float dR_lW_Z = std::sqrt(dPhi_lW_Z*dPhi_lW_Z + dEta_lW_Z*dEta_lW_Z);
        histMap1D_[{"dR_lW_Z", channel_, variation.first}]->Fill(dR_lW_Z, weight);
    }
    if (histMap1D_[{"nJetCSVv2T", channel_, variation.first}] != nullptr) {
        b_jetCSVv2->GetEntry(entry);
        unsigned int bjets = 0;
        for (const auto& jetCSVval : *jetCSVv2) {
            if (jetCSVval > 0.9535)
                bjets++;
        }
        histMap1D_[{"nJetCSVv2T", channel_, variation.first}]->Fill(bjets, weight);
    }
    if (histMap1D_[{"Pt", channel_, variation.first}] != nullptr) {
        b_Pt->GetEntry(entry);
        histMap1D_[{"Pt", channel_, variation.first}]->Fill(Pt, weight);
    }
    if (histMap1D_[{"nvtx", channel_, variation.first}] != nullptr) {
        b_nvtx->GetEntry(entry);
        histMap1D_[{"nvtx", channel_, variation.first}]->Fill(nvtx, weight);
    }
    if (histMap1D_[{"ZPt", channel_, variation.first}] != nullptr) {
        b_ZPt->GetEntry(entry);
        histMap1D_[{"ZPt", channel_, variation.first}]->Fill(ZPt, weight);
    }

    SafeHistFill(histMap1D_, "MTWZ", channel_, variation.first, MtWZ, weight*(isMC_ || MtWZ < 300 || noBlind));
    SafeHistFill(histMap1D_, "M3lMET", channel_, variation.first, M3lMET, weight*(isMC_ || M3lMET < 300 || noBlind));
    if (isMC_)
        SafeHistFill(histMap1D_, "nTruePU", channel_, variation.first, nTruePU, weight);
}

void WZSelector::SetupNewDirectory()
{
    SelectorBase::SetupNewDirectory();
    isaQGC_ = name_.find("aqgc") != std::string::npos;
    applyFullSelection_ = (selection_ == VBSselection_Loose_Full ||
                      selection_ == VBSselection_Tight_Full || 
                      selection_ == VBSselection_NoZeppenfeld_Full || 
                      selection_ == Inclusive2Jet_Full ||
                      selection_ == Wselection_Full ||
                      selection_ == VBSBackgroundControl_Full ||
                      selection_ == VBSBackgroundControlLoose_Full);
    //doSystematics_ = applyFullSelection_;
    //doSystematics_ = false;
    
    InitializeHistogramsFromConfig();   
}
