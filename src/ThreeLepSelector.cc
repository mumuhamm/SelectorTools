#include "Analysis/SelectorTools/interface/ThreeLepSelector.h"
#include "Analysis/SelectorTools/interface/HelperFunctions.h"

#include <TStyle.h>
#include <regex>
#include "TParameter.h"

#define Fill1D(NAME, VALUE_) HistFullFill(histMap1D_, NAME, variation.first, VALUE_, weight);
#define Fill2D(NAME, VALUE1_, VALUE2_) HistFullFill(histMap2D_, NAME, variation.first, VALUE1_, VALUE2_, weight);

#define GETMASK(index, size) (((1 << (size)) - 1) << (index))
#define READFROM(data, index, size) (((data) & GETMASK((index), (size))) >> (index))

#define CHGPT(index) (Electron_eCorr[index])
#define CLOSEJET_REWEIGHT

typedef std::bitset<sizeof(int)> IntBits;

enum ElectronCBID {CBID_VETO=1, CBID_LOOSE=2, CBID_MEDIUM=3, CBID_TIGHT=4};

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;

void ThreeLepSelector::SetScaleFactors() {
    calib = BTagCalibration("deepcsv", "data/btag_scales.csv");
    btag_reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
    btag_reader.load(calib, BTagEntry::FLAV_B, "comb");
    btag_reader.load(calib, BTagEntry::FLAV_C, "comb");
    btag_reader.load(calib, BTagEntry::FLAV_UDSG, "incl");

    pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
    if (pileupSF_ == nullptr )
	std::cout << "missing Pileup SF" << std::endl;

    // TFile* bFile = (TFile*) GetInputList()->FindObject("BScales");
    // name_ = GetInputList()->FindObject("name")->GetTitle();

    // Beff_j = (TH2D*) bFile->Get((name_ + "/Beff_j").c_str());
    // Beff_b = (TH2D*) bFile->Get((name_ + "/Beff_b").c_str());

    eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronTightIdSF");
    if (eIdSF_ == nullptr )
	std::cout  << "missing Electron ID SF" << std::endl;

    mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonMediumIdSF");
    if (mIdSF_ == nullptr )
	std::cout  << "missing Muon Id SF" << std::endl;

    TFile* f_btag_eff = new TFile("ScaleFactors/btageff__ttbar_powheg_pythia8_25ns_Moriond17_deepCSV.root");
    h_btag_eff_b    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_b");
    h_btag_eff_c    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_c");
    h_btag_eff_udsg = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_udsg");

}

void ThreeLepSelector::Init(TTree *tree) {
    b.SetTree(tree);

    allChannels_ = {{mm, "mm"}, {ee, "ee"}, {em, "em"}, {all, "all"}, {lll, "lll"}};

    hists1D_ = {
		"CutFlow",      "ZMass",       "ptl1",     "etal1",    "ptl2",     "etal2",        "SR",
		"bjetpt",       "jetpt",       "nbjet",    "njet",     "nleps",    "CRW_nbjet",    "CRW_njet",
		"CRZ_nbjet",    "CRZ_njet",    "CRZ_HT",   "CRZ_Met",  "Met",      "HT",           "weight",
		"CRW_HT",       "CRW_Met",     "CRZ_ptl3", "sphericity", "centrality",
    };
    hists2D_ = {"bJetvsJets",    "Beff_b_btag", "Beff_j_btag", "Beff_b", "Beff_j"};

    SelectorBase::Init(tree);
}

void ThreeLepSelector::SetBranchesNanoAOD() {
    //  NECESSARY!!!!
    b.CleanUp();
    
    b.SetBranch("HLT_DoubleMu8_Mass8_PFHT300", HLT_DoubleMu8_Mass8_PFHT300);
    b.SetBranch("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300);
    b.SetBranch("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300);
    b.SetBranch("HLT_AK8PFJet450", HLT_AK8PFJet450);
    b.SetBranch("HLT_PFJet450", HLT_PFJet450);

    b.SetBranch("nElectron",                  nElectron);
    b.SetBranch("Electron_pt",                Electron_pt);
    b.SetBranch("Electron_eta",               Electron_eta);
    b.SetBranch("Electron_phi",               Electron_phi);
    b.SetBranch("Electron_charge",            Electron_charge);
    b.SetBranch("Electron_mass",              Electron_mass);
    b.SetBranch("Electron_miniPFRelIso_all",  Electron_miniPFRelIso_all);
    b.SetBranch("Electron_dxy",               Electron_dxy);
    b.SetBranch("Electron_dz",                Electron_dz);
    b.SetBranch("Electron_sip3d",             Electron_sip3d);
    b.SetBranch("Electron_lostHits",          Electron_lostHits);
    b.SetBranch("Electron_convVeto",          Electron_convVeto);
    if(year_ == yr2018) {
	b.SetBranch("Electron_mvaFall17V2noIso",     Electron_MVA);
	b.SetBranch("Electron_cutBased",             Electron_cutBased);
    } else if(year_ == yr2017) {
	b.SetBranch("Electron_mvaFall17noIso",     Electron_MVA);
	b.SetBranch("Electron_cutBased",   Electron_cutBased);
    } else if(year_ == yr2016 || year_ == yrdefault) {
	b.SetBranch("Electron_mvaSpring16GP",        Electron_MVA);
	b.SetBranch("Electron_cutBased_Sum16",       Electron_cutBased);
    }
    b.SetBranch("Electron_tightCharge", Electron_tightCharge);
    b.SetBranch("Electron_sieie",                        Electron_sieie);
    b.SetBranch("Electron_hoe",                          Electron_hoe);
    b.SetBranch("Electron_deltaEtaSC",                   Electron_deltaEtaSC);
    b.SetBranch("Electron_eInvMinusPInv",                Electron_eInvMinusPInv);
    b.SetBranch("Electron_dr03EcalRecHitSumEt",          Electron_dr03EcalRecHitSumEt);
    b.SetBranch("Electron_dr03HcalDepth1TowerSumEt",     Electron_dr03HcalDepth1TowerSumEt);
    b.SetBranch("Electron_dr03TkSumPt",                  Electron_dr03TkSumPt);
    b.SetBranch("Electron_vidNestedWPBitmapSum16", Electron_vidBitmap);
    b.SetBranch("Electron_eCorr", Electron_eCorr);
    
    b.SetBranch("nMuon",                  nMuon);
    b.SetBranch("Muon_pt",                Muon_pt);
    b.SetBranch("Muon_eta",               Muon_eta);
    b.SetBranch("Muon_phi",               Muon_phi);
    b.SetBranch("Muon_mass",              Muon_mass);
    b.SetBranch("Muon_charge",            Muon_charge);
    b.SetBranch("Muon_mediumId",          Muon_mediumId);
    b.SetBranch("Muon_miniPFRelIso_all",  Muon_miniPFRelIso_all);
    b.SetBranch("Muon_dxy",               Muon_dxy);
    b.SetBranch("Muon_dz",                Muon_dz);
    b.SetBranch("Muon_sip3d",             Muon_sip3d);
    b.SetBranch("Muon_isGlobal",          Muon_isGlobal);
    b.SetBranch("Muon_isTracker",          Muon_isTracker);
    b.SetBranch("Muon_isPFcand",          Muon_isPFcand);
    b.SetBranch("Muon_tightCharge",       Muon_tightCharge);
    
    b.SetBranch("nJet",               nJet);
    b.SetBranch("Jet_btagCSVV2",      Jet_btagCSVV2);
    b.SetBranch("Jet_btagDeepB",      Jet_btagDeepB);
    b.SetBranch("Jet_eta",            Jet_eta);
    b.SetBranch("Jet_phi",            Jet_phi);
    b.SetBranch("Jet_pt",             Jet_pt);
    b.SetBranch("Jet_mass",           Jet_mass);
    b.SetBranch("Jet_jetId",          Jet_jetId);
    b.SetBranch("Jet_hadronFlavour",  Jet_hadronFlavour);
    b.SetBranch("Jet_rawFactor",  Jet_rawFactor);

    b.SetBranch("MET_pt",     MET);
    b.SetBranch("MET_phi",    type1_pfMETPhi);

    b.SetBranch("event", event);
    b.SetBranch("luminosityBlock", lumi);

    if(applyScaleFactors_) {
	b.SetBranch("GenPart_pdgId",GenPart_pdgId);
	b.SetBranch("GenPart_genPartIdxMother", GenPart_genPartIdxMother);
	b.SetBranch("GenPart_status", GenPart_status);
	b.SetBranch("nGenPart", nGenPart);
    }
    
    b.SetBranch("Flag_goodVertices", Flag_goodVertices);
    b.SetBranch("Flag_globalSuperTightHalo2016Filter", Flag_globalSuperTightHalo2016Filter);
    b.SetBranch("Flag_HBHENoiseFilter", Flag_HBHENoiseFilter);
    b.SetBranch("Flag_HBHENoiseIsoFilter", Flag_HBHENoiseIsoFilter);
    b.SetBranch("Flag_EcalDeadCellTriggerPrimitiveFilter", Flag_EcalDeadCellTriggerPrimitiveFilter);
    b.SetBranch("Flag_BadPFMuonFilter", Flag_BadPFMuonFilter);
    b.SetBranch("Flag_ecalBadCalibFilter", Flag_ecalBadCalibFilter);

#ifdef CLOSEJET_REWEIGHT
    b.SetBranch("Jet_L1", Jet_L1);
    b.SetBranch("Jet_L2L3", Jet_L2L3);
#endif // CLOSEJET_REWEIGHT

    if (isMC_) {
	b.SetBranch("genWeight",    genWeight);
	b.SetBranch("Pileup_nPU",   numPU);
	b.SetBranch("Pileup_nTrueInt", Pileup_nTrueInt);
    }

    mvaValues[CBID_LOOSE] = {{-0.46, -0.48, 0, -0.85}, 
			     {-0.03, -0.67, 0, -0.91},
			     {0.06, -0.49, 0, -0.83}};
    mvaValues[CBID_TIGHT] = {{10, 0.77, 0, 0.52}, 
			     {10, 0.56, 0, 0.11},
			     {10, 0.48, 0, -0.01}};
    /// Setup interpolation used for 25>pt>15
    for(auto& pair: mvaValues) {
	for(auto& vec: pair.second) {
	    vec[2] = (vec[3]-vec[1])/10;
	}
    }
}

/// Make to seperate fuctionality
void ThreeLepSelector::clearValues() {
    weight = 1;
    HT = 0;
    nJets = 0;
    nBJets = 0;
    passZVeto = true;
    goodLeptons.clear();
    looseLeptons.clear();
    goodJets.clear();
}

void ThreeLepSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) {
    clearValues();

    b.SetEntry(entry);

    if (nElectron > N_KEEP_MU_E_ || nMuon > N_KEEP_MU_E_ || nGenPart > N_KEEP_GEN_) {
	std::string message = "Found more electrons or muons than max read number.\n    Found ";
	message += std::to_string(nElectron);
	message += " electrons.\n    Found ";
	message += std::to_string(nMuon);
	message += " Muons\n  --> Max read number was ";
	message += std::to_string(N_KEEP_MU_E_);
	message += "\nExiting because this can cause problems. Increase N_KEEP_MU_E_ to avoid this error.\n";
	message += std::to_string(nGenPart);
	message += " Gens\n  --> Max read number was ";
	message += std::to_string(N_KEEP_GEN_);
	message += "\nExiting because this can cause problems. Increase N_KEEP_GEN_ to avoid this error.\n";
	throw std::domain_error(message);
    }

    if(eventVec[event]) return;
    
    /// basic setups
    setupElectrons();
    setupMuons();
    setupJets();
    setupChannel();

    if (isMC_) {
	ApplyScaleFactors();
    }

    if(goodLeptons.size() >= 2 && looseLeptons.size() >= 3) {
	for(auto lep : goodLeptons) {
	    passZVeto = doesPassZVeto(lep, looseLeptons);
	    if(!passZVeto) break;
	}
    }
}

void ThreeLepSelector::setupMuons() {
    for (size_t i = 0; i < nMuon; ++i) {
	if(isGoodMuon(i)) {
	    goodLeptons.push_back(GoodPart(get4Vector(PID_MUON, i), PID_MUON * Muon_charge[i]));
	    goodLeptons.back().index = i;
	    looseLeptons.push_back(goodLeptons.back());
            
	    if(!passFullIso(goodLeptons.back().v, 0.76, 7.2)) {    // Extra Iso requirement
		goodLeptons.pop_back();
	    }
	}
	else if(isLooseMuon(i)) {
	    looseLeptons.push_back(GoodPart(get4Vector(PID_MUON, i), PID_MUON * Muon_charge[i]));
	    looseLeptons.back().index = i;
	}
    }

}

void ThreeLepSelector::setupElectrons() {
    for (size_t i = 0; i < nElectron; ++i) {
	if( isGoodElectron(i)) {
            goodLeptons.push_back(GoodPart(get4Vector(PID_ELECTRON, i), PID_ELECTRON * Electron_charge[i]));
	    goodLeptons.back().index = i;
	    looseLeptons.push_back(goodLeptons.back());
            
	    if(!passFullIso(goodLeptons.back().v, 0.8, 7.2)) {   // Extra Iso requirement
	    	goodLeptons.pop_back();
	    }
	}
	else if(isLooseElectron(i)) {
	    looseLeptons.push_back(GoodPart(get4Vector(PID_ELECTRON, i), PID_ELECTRON * Electron_charge[i]));
	    looseLeptons.back().index = i;
	}
    }
}

void ThreeLepSelector::setupJets() {
    std::vector<size_t> closeJet;
    for(auto lep : looseLeptons) {
    	if(passFakeableCuts(lep)) closeJet.push_back(getCloseJetIndex(lep.v, 0.4));
    }

    for(size_t i = 0; i < nJet; ++i) {
	// if(goodLeptons.size() < 2) break;  // only try to find jets if have leptons
	if(std::find(closeJet.begin(), closeJet.end(), i) != closeJet.end()) continue;
	/// jet
	bool passedGoodJet = false;
	if(isGoodJet(i)) {
	    passedGoodJet = true;
	    goodJets.push_back(GoodPart(get4Vector(PID_JET, i), Jet_hadronFlavour[i]));
	    goodJets.back().setPassJetSel();
	    goodJets.back().index = i;
	    nJets++;
	    HT += Jet_pt[i];
	}
	// bjet
	if(isGoodBJet(i)) {
	    if(!passedGoodJet) {
		goodJets.push_back(GoodPart(get4Vector(PID_BJET, i), Jet_hadronFlavour[i]));
		goodJets.back().index = i;
	    }
            goodJets.back().setPassBJetSel();
	    nBJets++;
	}
    }

}

void ThreeLepSelector::setupChannel() {
    if(goodLeptons.size() > 2) {
	channelName_ = "lll";
	if(goodLeptons[1].Charge() * goodLeptons[2].Charge() > 0) {
	    std::swap(goodLeptons[0], goodLeptons[2]);
	}
	else if(goodLeptons[0].Charge() * goodLeptons[2].Charge() > 0) {
	    std::swap(goodLeptons[1], goodLeptons[2]);
	}
	/// PT swap
	if(goodLeptons[0].Pt() < goodLeptons[1].Pt()) {
	    std::swap(goodLeptons[0], goodLeptons[1]);
	}
    }
    else if(goodLeptons.size() < 2) {
	channelName_ = "Unknown";
    }
    else if(goodLeptons[0].Id() == PID_MUON && goodLeptons[1].Id() == PID_MUON) {
	channelName_ = "mm";
    }
    else if(goodLeptons[0].Id() == PID_ELECTRON && goodLeptons[1].Id() == PID_ELECTRON) {
	channelName_ = "ee";
    }
    else {
	channelName_ = "em";
	if(goodLeptons[0].Pt() < goodLeptons[1].Pt()) {
	    std::swap(goodLeptons[0], goodLeptons[1]);
	}
    }

    channel_ = channelMap_[channelName_];
}


LorentzVector ThreeLepSelector::get4Vector(PID pid, int idx) {
    if(pid == PID_MUON)
	return LorentzVector(Muon_pt[idx], Muon_eta[idx], Muon_phi[idx], Muon_mass[idx]);
    else if(pid == PID_ELECTRON)
	return LorentzVector(Electron_pt[idx]/CHGPT(idx), Electron_eta[idx], Electron_phi[idx], Electron_mass[idx]);
    else
	return LorentzVector(Jet_pt[idx], Jet_eta[idx], Jet_phi[idx], Jet_mass[idx]);
}


bool ThreeLepSelector::doesPassZVeto(GoodPart& lep, std::vector<GoodPart>& looseList) {
    for (auto lLep : looseList) {
	if( lep.Charge() == -1*lLep.Charge() &&  //opposite charge, same ID
	    (abs((lLep.v + lep.v).M() - 91.188) < 15 || (lLep.v + lep.v).M() < 12)) {
	    return false;
	}
    }
    return true;
}

void ThreeLepSelector::ApplyScaleFactors() {
    if(selection_ == BEfficiency) return;
    // weight *= (genWeight > 0) ? 1 : -1;
    weight *= genWeight;

    if(!applyScaleFactors_ || goodLeptons.size() < 2) return;

    for(auto lep: goodLeptons) {
	weight *= leptonScaleFactor(lep.Id(), lep.Pt(), lep.Eta(), HT);
    }
        
    weight *= triggerScaleFactor(goodLeptons[0].Id(), goodLeptons[1].Id(),
				 goodLeptons[0].Pt(), goodLeptons[1].Pt(),
				 goodLeptons[0].Eta(), goodLeptons[1].Eta(), HT);
    weight *= getTruePUw_Moriond(Pileup_nTrueInt);
    weight *= getWDecayScaleFactor();
    
    for(auto jet : goodJets) {
	BTagEntry::JetFlavor flav;
	if(jet.Id() == PID_BJET) flav = BTagEntry::FLAV_B;
	else if(jet.Id() == PID_CJET) flav = BTagEntry::FLAV_C;
	else  flav = BTagEntry::FLAV_UDSG;
	double bSF = btag_reader.eval_auto_bounds("central",  flav, jet.Eta(), jet.Pt());
	if( jet.passedBJetSel() ) {
	    weight *= bSF;
	} else {
	    double eff = getBtagEffFromFile(jet.Pt(), jet.Eta(), jet.Id());
	    weight *= (1 - bSF * eff) / (1 - eff);
	}
    }

    return;
}
bool ThreeLepSelector::passFakeableCuts(GoodPart& lep) {
    int index = lep.index;
    if(lep.Pt() < 10) return false;
    if(lep.Id() == PID_MUON) {
	return (Muon_mediumId[index] 
		&& Muon_sip3d[index] < 4
		&& Muon_tightCharge[index] == 2
		&& passFullIso(lep.v, 0.72, 7.2)
		
		);
    }
    else {
	return (Electron_sip3d[index] < 4
		&& Electron_tightCharge[index] == 2
		&& Electron_lostHits[index] == 0
		&& passFullIso(lep.v, 0.8, 7.2)
		);
    }
}


bool ThreeLepSelector::isGoodMuon(size_t index) {
    bool yearCuts = true;
    if(year_ == yr2016) yearCuts = (Muon_miniPFRelIso_all[index] < 0.16);
    else                yearCuts = (Muon_miniPFRelIso_all[index] < 0.11);
    
    return ( (Muon_pt[index] > 20) 
	     && (Muon_tightCharge[index] == 2) 
	     && (abs(Muon_eta[index]) < 2.4) 
	     && (Muon_mediumId[index]) 
	     && (yearCuts) 
	     && (abs(Muon_dz[index]) < 0.1) 
	     && (abs(Muon_dxy[index]) < 0.05) 
	     && (Muon_sip3d[index] < 4)
	     );
}

bool ThreeLepSelector::passMVACut(std::vector<std::vector<double> > mvaCuts, int index) {
    int caseIndex = 0;
    //// PT Splitting
    if(Electron_pt[index]/CHGPT(index) < 5)       return false;
    else if(Electron_pt[index]/CHGPT(index) < 10) caseIndex += 0;
    else if(Electron_pt[index]/CHGPT(index) < 15) caseIndex += 1;
    else if(Electron_pt[index]/CHGPT(index) < 25) caseIndex += 2;
    else                             caseIndex += 3;
    //// ETA Splitting
    if(abs(Electron_eta[index]) < 0.8)        caseIndex += 0;
    else if(abs(Electron_eta[index]) < 1.479) caseIndex += 4;
    else if(abs(Electron_eta[index]) < 2.5)   caseIndex += 8;

    if(caseIndex % 4 != 2) return Electron_MVA[index] > mvaCuts[caseIndex/4][caseIndex%4];
    else                  return Electron_MVA[index] > mvaInterpolate(Electron_pt[index]/CHGPT(index), mvaCuts[caseIndex/4]);
}

double ThreeLepSelector::mvaInterpolate(double pt, std::vector<double> cuts) {
    return cuts[1] + cuts[2]*(pt-15);
}


bool ThreeLepSelector::isGoodElectron(size_t index) {
    if(abs(Electron_eta[index]) > 2.5) return false;
    bool passId = false;
    
    if(selection_ == FourTopMVAEl || selection_ != FourTopCutBasedEl) {
	if(year_ == yr2016 || year_ == yrdefault) {
            passId = passMVACut(mvaValues[CBID_TIGHT], index);
	    passId = passId && (Electron_miniPFRelIso_all[index] < 0.12);
	}
	else if(year_ == yr2017) {
	    ///// NEED to fix mva values for 2017
	    passId = passMVACut(mvaValues[CBID_TIGHT], index);
	    passId = passId && (Electron_miniPFRelIso_all[index] < 0.07);
	}
    } else {
	passId = (Electron_cutBased[index] >= CBID_MEDIUM);
	//	if(year_ == yr2016)       passId = passId && (Electron_miniPFRelIso_all[index] < 0.12);
	//	else if(year_ == yr2017)       passId = passId && (Electron_miniPFRelIso_all[index] < 0.07);
    }

    return ((Electron_pt[index]/CHGPT(index) > 20)
	    && (passId)
	    && (Electron_convVeto[index]) 
	    && (Electron_lostHits[index] == 0) 
	    && (Electron_tightCharge[index] == 2) 
	    && (abs(Electron_dz[index]) < 0.1) 
	    && (abs(Electron_dxy[index]) < 0.05) 
	    && (Electron_sip3d[index] < 4)
	    && passTriggerEmu(index)
	    && Electron_dr03EcalRecHitSumEt[index] / Electron_pt[index]*CHGPT(index) < 0.45 
	    && Electron_dr03HcalDepth1TowerSumEt[index] / Electron_pt[index]*CHGPT(index) < 0.25 
	    && Electron_dr03TkSumPt[index] / Electron_pt[index]*CHGPT(index) < 0.2 
	    );
}


bool ThreeLepSelector::isLooseMuon(size_t index) {
    return ((Muon_isGlobal[index] || Muon_isTracker[index])
	    && (Muon_isPFcand[index])
	    && (Muon_miniPFRelIso_all[index] < 0.4)
	    && (abs(Muon_dz[index]) < 0.1)
	    && (abs(Muon_dxy[index]) < 0.05)
	    && (Muon_pt[index] > 5)
	    );
}

bool ThreeLepSelector::isLooseElectron(size_t index) {
    bool passId = false;

    if(selection_ == FourTopMVAEl || selection_ != FourTopCutBasedEl) {
	passId = passMVACut(mvaValues[CBID_LOOSE], index);
    }
    else {
	passId = (Electron_cutBased[index] >= CBID_LOOSE);
    }
    return ((passId)//
	    && (Electron_pt[index]/CHGPT(index) > 7)
	    && (Electron_convVeto[index]) 
	    && (Electron_lostHits[index] <= 1)
	    && (Electron_miniPFRelIso_all[index] < 0.4) 
	    && passTriggerEmu(index)
	    && (abs(Electron_dz[index]) < 0.1) 
	    && (abs(Electron_dxy[index]) < 0.05)
	    );
}

bool ThreeLepSelector::isGoodJet(size_t index) {
    bool yearCut = true;
    double ptCut = 40;
    
    if(year_ == yr2016) yearCut = IntBits(Jet_jetId[index]).test(0);
    else                yearCut = IntBits(Jet_jetId[index]).test(1);
    return ((Jet_pt[index] > ptCut)      &&
	    (abs(Jet_eta[index]) < 2.4) &&
	    (yearCut)
	    );
}

/// TODO: add toggle for different btag stuff
// (Jet_btagCSVV2[index] > 0.8484) &&
bool ThreeLepSelector::isGoodBJet(size_t index) {
    bool yearCut = true;
    double ptCut = 25;
    
    if(year_ == yr2016) yearCut = (IntBits(Jet_jetId[index]).test(0)) && (Jet_btagDeepB[index] > 0.6324);
    else                yearCut = (IntBits(Jet_jetId[index]).test(1)) && (Jet_btagDeepB[index] > 0.4941);

    return ((Jet_pt[index] > ptCut)
	    && (abs(Jet_eta[index]) < 2.4)
	    && (yearCut)
	    );
}

size_t ThreeLepSelector::getCloseJetIndex(LorentzVector& lep, double minDR ) {
    size_t minIndex = -1;
    
    for(size_t index = 0; index < nJet; ++index) {
	LorentzVector jet = get4Vector(PID_JET, index);
	double dr = reco::deltaR(jet, lep);
	if(minDR > dr) {
	    minDR = dr;
	    minIndex = index;
	}
    }
    return minIndex;
}

bool ThreeLepSelector::passFullIso(LorentzVector& lep, double I2, double I3) {
    int closeIdx = getCloseJetIndex(lep);
    LorentzVector closeJet  = get4Vector(PID_JET, closeIdx);
#ifdef CLOSEJET_REWEIGHT
    closeJet = (Jet_L1[closeIdx]*(1-Jet_rawFactor[closeIdx])*closeJet-lep)*Jet_L2L3[closeIdx]+lep;
#endif // CLOSEJET_REWEIGHT

    return (lep.Pt()/closeJet.Pt() > I2) || (LepRelPt(lep, closeJet) > I3);
}


double ThreeLepSelector::LepRelPt(LorentzVector& lep, LorentzVector& closeJet) {
    auto diff = closeJet.Vect() - lep.Vect();
    auto cross = diff.Cross(lep.Vect());
    return std::sqrt(cross.Mag2()/diff.Mag2());
}

// Need to include DeltaPhi Ieta
bool ThreeLepSelector::passTriggerEmu(size_t index) {
    bool etaDepend = false;
    int DEtaInCut = READFROM(Electron_vidBitmap[index], 6, 3);
    int DPhiInCut = READFROM(Electron_vidBitmap[index], 9, 3);

    if(abs(Electron_eta[index]) < 1.479) {
	etaDepend = Electron_sieie[index] < 0.011 && (DEtaInCut >= 1) && (DPhiInCut >= 4);
    }
    else {
	etaDepend = Electron_sieie[index] < 0.031 && (DEtaInCut >= 1) && (DPhiInCut >= 3);
    }
    return (abs(Electron_eInvMinusPInv[index]) < 0.01 &&
	    etaDepend &&
	    Electron_hoe[index] < 0.08
	    );
}


bool ThreeLepSelector::doesNotOverlap(size_t index) {
    LorentzVector tmp = get4Vector(PID_JET, index);
    double dR = 0.4;
    for(auto lep: goodLeptons) {
	if(reco::deltaR(tmp, lep.v) < dR) return false;
    }
    return true;
}

bool ThreeLepSelector::MetFilter() {
    return Flag_goodVertices
	&& Flag_globalSuperTightHalo2016Filter
	&& Flag_HBHENoiseFilter
	&& Flag_HBHENoiseIsoFilter
	&& Flag_EcalDeadCellTriggerPrimitiveFilter
	&& Flag_BadPFMuonFilter
	&& Flag_ecalBadCalibFilter;
}

void ThreeLepSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) {
    int step = 0;
    Fill1D("CutFlow", 0);

    /// Trigger
    if(((channel_ == mm && !HLT_DoubleMu8_Mass8_PFHT300) ||
	(channel_ == em && !HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300) ||
	(channel_ == ee && !HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300))
       && !HLT_AK8PFJet450 && !HLT_PFJet450
       ) return;
    Fill1D("CutFlow", ++step);

    if(!MetFilter()) return;
    Fill1D("CutFlow", ++step);

    /// 2 good leptons
    if(goodLeptons.size() < 2) return;
    Fill1D("CutFlow", ++step);

    // first lep requirement
    if(goodLeptons[0].Pt() < 25) return;
    Fill1D("CutFlow", ++step);

    // same sign requirement
    if((goodLeptons.size() == 2 && goodLeptons[0].Charge() * goodLeptons[1].Charge() < 0)
       || (goodLeptons.size() == 3 && goodLeptons[0].Charge() * goodLeptons[2].Charge() > 0)
       )
	return;
    Fill1D("CutFlow", ++step);

    // jet cut
    if(nJets < 2) return;
    Fill1D("CutFlow", ++step);

    // bjet cut
    if(nBJets < 2) return;
    Fill1D("CutFlow", ++step);

    // ht cut
    if(HT < 300) return;
    Fill1D("CutFlow", ++step);

    // met cut
    if (MET < 50) return;
    Fill1D("CutFlow", ++step);

    Fill1D("SR", getSRBin());

    if(getSRBin() == -1) {
	return;
    }
    else if(getSRBin() == 0) {
	Fill1D("CRW_Met", MET);
	Fill1D("CRW_HT", HT);
	Fill1D("CRW_njet", nJets);
	Fill1D("CRW_nbjet", nBJets);
	return;
    }
    else if(getSRBin() == 9) {
	Fill1D("CRZ_Met", MET);
	Fill1D("CRZ_HT", HT);
    	Fill1D("CRZ_njet", nJets);
	Fill1D("CRZ_nbjet", nBJets);
	Fill1D("CRZ_ptl3", goodLeptons[3].Pt());
	return;
    }

    Fill1D("CutFlow", ++step);
    
    HistFullFill(histMap1D_, "weight", variation.first, abs(weight), 1);
    Fill1D("Met", MET);
    Fill1D("HT", HT);
    Fill1D("ptl1", goodLeptons[0].Pt());
    Fill1D("ptl2", goodLeptons[1].Pt());
    Fill1D("njet", nJets);
    Fill1D("nbjet", nBJets);
    Fill2D("bJetvsJets", nJets, nBJets);
    Fill1D("nleps", goodLeptons.size());

    Fill1D("sphericity", JetSphericity(goodJets));
    Fill1D("centrality", JetCentrality(goodJets,HT));
    
    
}

void ThreeLepSelector::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    InitializeHistogramsFromConfig();
}

int ThreeLepSelector::getSRBin() const {
    if(goodLeptons.size() == 2) {
	if(nBJets == 2) {
	    if(!passZVeto)      return -1;
	    else if(nJets <= 5)   return 0;  // WCR
	    else if(nJets == 6)  return 1;
	    else if(nJets == 7)  return 2;
	    else if(nJets >= 8)   return 3;
	}
	else if(nBJets == 3) {
	    if(nJets == 5)       return 4;
	    else if(nJets == 6)  return 4;
	    else if(nJets >= 7)   return 5;
	}
	else if(nBJets >= 4) {
	    if(nJets >= 5)         return 6;
	}
    } else {
	if(!passZVeto)                     return 9;  /// ZCR
	else if(nBJets == 2 && nJets >= 5)   return 7;
	else if(nBJets >= 3 && nJets >= 4)    return 8;
    }

    return -1;
}

float ThreeLepSelector::getBtagEffFromFile(double pt, double eta, int mcFlavour){
    float pt_cutoff = std::max(20.,std::min(399.,pt));
    if (abs(mcFlavour) == 5) {
	// use pt bins up to 600 GeV for b
	pt_cutoff = std::max(20.,std::min(599.,pt));
	return h_btag_eff_b->GetBinContent(h_btag_eff_b->FindBin(pt_cutoff, fabs(eta)));
    }
    else if (abs(mcFlavour) == 4) {
	return h_btag_eff_c->GetBinContent(h_btag_eff_c->FindBin(pt_cutoff, fabs(eta)));
    }
    else {
	return h_btag_eff_udsg->GetBinContent(h_btag_eff_udsg->FindBin(pt_cutoff, fabs(eta)));
    }
}

double ThreeLepSelector::getWDecayScaleFactor() {
    float pdgleptW = 0.3258;
    float genleptW = 1.0/3;
    int nleptonicW = 0;
    int nW = 0;
    
    for(size_t i=0; i < nGenPart; i++) {
    	if(abs(GenPart_pdgId[i]) == 24
    	   && (GenPart_status[i] == 22 || GenPart_status[i] == 52)
    	   && abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) != 24) {
    	    nW++;
    	}
    	else if((abs(GenPart_pdgId[i]) == 12
    		 || abs(GenPart_pdgId[i]) == 14
    		 || abs(GenPart_pdgId[i]) == 16)
    		&& abs(GenPart_pdgId[GenPart_genPartIdxMother[i]]) == 24) {
    	    nleptonicW++;
    	}
    	else continue;
    }

    int nhadronicW = nW - nleptonicW;
    return pow((pdgleptW/genleptW),nleptonicW) * pow(((1-pdgleptW)/(1-genleptW)),nhadronicW);
}
