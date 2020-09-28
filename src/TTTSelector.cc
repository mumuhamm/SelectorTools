#include "Analysis/SelectorTools/interface/TTTSelector.h"

#include <TStyle.h>
#include <regex>
#include "TParameter.h"

#define Fill1D(NAME, VALUE_) HistFullFill(histMap1D_, NAME, variation.first, VALUE_, weight);
#define Fill2D(NAME, VALUE1_, VALUE2_) HistFullFill(histMap2D_, NAME, variation.first, VALUE1_, VALUE2_, weight);



enum Lepton {Muon=13, Electron=11};

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LorentzVector;

void TTTSelector::SlaveBegin(TTree *) {
  return;
  pileupSF_ = (ScaleFactor *) GetInputList()->FindObject("pileupSF");
  if (pileupSF_ == nullptr ) 
    Abort("Must pass pileup weights SF");
  eIdSF_ = (ScaleFactor *) GetInputList()->FindObject("electronTightIdSF");
  if (eIdSF_ == nullptr ) 
    Abort("Must pass electron ID SF");
  eGsfSF_ = (ScaleFactor *) GetInputList()->FindObject("electronGsfSF");
  if (eGsfSF_ == nullptr ) 
    Abort("Must pass electron GSF SF");
  mIdSF_ = (ScaleFactor *) GetInputList()->FindObject("muonTightIdSF");
  if (mIdSF_ == nullptr ) 
    Abort("Must pass muon ID SF");
  mIsoSF_ = (ScaleFactor *) GetInputList()->FindObject("muonIsoSF");
  if (mIsoSF_ == nullptr ) 
    Abort("Must pass muon Iso SF");

  prefireEff_ = (TEfficiency*) GetInputList()->FindObject("prefireEfficiencyMap");
  if (prefireEff_ == nullptr ) 
    Abort("Must pass prefiring efficiency map");

}

void TTTSelector::Init(TTree *tree){
  b.SetTree(tree);
  
  allChannels_ = {{ee, "ee"}, {mm, "mm"}, {em, "em"}, {all, "all"}};
 
  hists1D_ = {"CutFlow", "ZMass", "ptl1", "etal1", "ptl2", "etal2", "SR", "bjetpt", "jetpt", "nbjet", "njet", "jetphi","bjetphi", "phil1", "phil2","HT","MET", "nelec", "nmuon", "lept_charge", "before", "after", "nGoodParts"};
 hists2D_ = {"bJetvsJets"};
  
  SelectorBase::Init(tree);
  
}

void TTTSelector::SetBranchesNanoAOD() {
  //  NECESSARY!!!!
  b.CleanUp();

  b.SetBranch("nElectron", nElectron);
  b.SetBranch("Electron_pt", Electron_pt);
  b.SetBranch("Electron_eta", Electron_eta);
  b.SetBranch("Electron_phi", Electron_phi);
  b.SetBranch("Electron_charge", Electron_charge);
  b.SetBranch("Electron_mass", Electron_mass);
  b.SetBranch("Electron_mass", Electron_miniPFRelIso_all);
  b.SetBranch("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all);
  b.SetBranch("Electron_dxy", Electron_dxy);
  b.SetBranch("Electron_dz", Electron_dz);
  b.SetBranch("Electron_sip3d", Electron_sip3d);
  if(year_ == yr2018) {
    b.SetBranch("Electron_mvaFall17V2noIso", Electron_MVA);
    b.SetBranch("Electron_cutBased", Electron_cutBased);
  } else if(year_ == yr2017) {
    b.SetBranch("Electron_mvaFall17V1noIso", Electron_MVA);
    b.SetBranch("Electron_cutBased_Fall17_V1", Electron_cutBased);
  } else if(year_ == yr2016 || year_ == yrdefault) {
    b.SetBranch("Electron_mvaSpring16GP", Electron_MVA);
    b.SetBranch("Electron_cutBased_Sum16", Electron_cutBased);
  }

  b.SetBranch("nMuon", nMuon);
  b.SetBranch("Muon_pt", Muon_pt);
  b.SetBranch("Muon_eta", Muon_eta);
  b.SetBranch("Muon_phi", Muon_phi);
  b.SetBranch("Muon_tightId", Muon_tightId);
  b.SetBranch("Muon_mediumId", Muon_mediumId);
  b.SetBranch("Muon_tkIsoId", Muon_tkIsoId);
  b.SetBranch("Muon_pfRelIso04_all", Muon_pfRelIso04_all);
  b.SetBranch("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all);
  b.SetBranch("Muon_charge", Muon_charge);
  b.SetBranch("Muon_mass", Muon_mass);
  b.SetBranch("Muon_dxy", Electron_dxy);
  b.SetBranch("Muon_dz", Electron_dz);
  b.SetBranch("Muon_sip3d", Electron_sip3d);
  
  b.SetBranch("nJet", nJet);
  b.SetBranch("Jet_btagCSVV2", Jet_btagCSVV2);
  b.SetBranch("Jet_btagDeepB", Jet_btagDeepB);
  b.SetBranch("Jet_eta", Jet_eta);
  b.SetBranch("Jet_phi", Jet_phi);
  b.SetBranch("Jet_pt", Jet_pt);
  b.SetBranch("Jet_mass", Jet_mass);

  b.SetBranch("Jet_neHEF", Jet_neHEF);
  b.SetBranch("Jet_neEmEF", Jet_neEmEF);
  b.SetBranch("Jet_nConstituents", Jet_nConstituents);
  b.SetBranch("Jet_chHEF", Jet_chHEF);
  b.SetBranch("Jet_chEmEF", Jet_chEmEF);

  b.SetBranch("MET_pt", MET);
  b.SetBranch("MET_phi", type1_pfMETPhi);
  
  if (isMC_) {
    b.SetBranch("genWeight", genWeight);
    b.SetBranch("Pileup_nPU", numPU);
  }
  
}

void TTTSelector::SetBranchesUWVV() {
  return;
}

void TTTSelector::LoadBranchesUWVV(Long64_t entry, std::pair<Systematic, std::string> variation) {
  return;
}

/// Make to seperate fuctionality
void TTTSelector::clearValues() {
  weight = 1;
  HT = 0;
  nTightJet = 0;
  nBJets = 0;
  goodParts.clear();
  charge.clear();
}

void TTTSelector::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
  clearValues();
  b.SetEntry(entry);

  if (nElectron > N_KEEP_MU_E_ || nMuon > N_KEEP_MU_E_) {
    std::string message = "Found more electrons or muons than max read number.\n    Found ";
    message += std::to_string(nElectron);
    message += " electrons.\n    Found ";
    message += std::to_string(nMuon);
    message += " Muons\n  --> Max read number was ";
    message += std::to_string(N_KEEP_MU_E_);
    message += "\nExiting because this can cause problems. Increase N_KEEP_MU_E_ to avoid this error.\n";
    throw std::domain_error(message);
  }

  /////////////////////
  // Setup Electrons //
  /////////////////////

  for (size_t i = 0; i < nElectron; i++) {
    // if(IsGoodMVAElectron2016(i)) {
    if(IsGoodCBElectron(i)) {
      // // Extra Iso requirement
      // LorentzVector lep(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
      // if(!passFullIso(lep, 0.8, 7.2)) continue;

      // Setup goodPart
      goodParts.push_back(GoodPart(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]));
      goodParts.back().SetPdgId(Electron*Electron_charge[i]);
      charge.push_back(Electron_charge[i]);
    }
  }

  /////////////////
  // Setup Muons //
  /////////////////
  
  for (size_t i = 0; i < nMuon; i++) {
    if(IsGoodMuon(i)) {
      // // Extra Iso requirement
      // LorentzVector lep(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
      // if(!passFullIso(lep, 0.76, 7.2)) continue;

      // Setup goodPart
      goodParts.push_back(GoodPart(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]));
      goodParts.back().SetPdgId(Muon*Muon_charge[i]);
      charge.push_back(Muon_charge[i]);
      
    }
  }
  
  ////////////////
  // Setup Jets //
  ////////////////
  
  for(size_t i = 0; i < nJet; i++) {
    if(goodParts.size() != 2) break;

    //common to both jets, may put more here
    if(!isLooseJetId(i)) continue;

    // regular jets
    if(IsGoodJet(i)) {
      nTightJet++;
      HT += Jet_pt[i];
    }

    // bjet 
    if(IsGoodBJet(i)) nBJets++;
  }

  channel_ = channelMap_[channelName_];
  if(goodParts.size() != 2) {
    channel_ = Unknown;
    channelName_ = "Unknown";
  } else if(goodParts[0].Id() == Muon && goodParts[1].Id() == Muon) {
    channel_ = mm;
    channelName_ = "mm";
  } else if(goodParts[0].Id() == Electron && goodParts[1].Id() == Electron) {
    channel_ = ee;
    channelName_ = "ee";
  } else {
    channel_ = em;
    channelName_ ="em";
    /// fix order of leptons by pt
    if(goodParts[0].Pt() < goodParts[1].Pt()) {
      std::swap(goodParts[0], goodParts[1]);
    }
  }

  if (isMC_) {
    ApplyScaleFactors();
  }
}


void TTTSelector::ApplyScaleFactors() {
  weight *= (genWeight > 0) ? 1.0 : -1.0;
  //  weight = 1;
  // skipping for now
  return;
  //// to add!!
}

////////////////////////////////////////////////////
// Functions for defining particles: used in main //
// loop of leptons/jets for getting multiplicity  //
////////////////////////////////////////////////////



bool TTTSelector::IsGoodMuon(size_t index) {
  return ( (Muon_pt[index] > 20) &&
	   (abs(Muon_eta[index]) < 2.4) &&
	   Muon_mediumId[index] &&
	   (Muon_miniPFRelIso_all[index] < 0.16) &&
	   (Muon_dz[index] < 0.1) &&
	   (Muon_dxy[index] < 0.05) &&
	   (Muon_sip3d[index] < 4));
}

bool TTTSelector::IsGoodCBElectron(size_t index) {
  return ((Electron_pt[index] > 20) &&
	  (abs(Electron_eta[index]) < 2.5) &&
	  (Electron_miniPFRelIso_all[index] < 0.12) &&
	  (Electron_cutBased[index] == 4)  &&
	  (Electron_dz[index] < 0.1) &&
	  (Electron_dxy[index] < 0.05) &&
	  (Electron_sip3d[index] < 4));
}

bool TTTSelector::IsGoodMVAElectron2016(size_t index) {
  bool mvaRec = false;
  if(abs(Electron_eta[index]) < 0.8)
    mvaRec = Electron_MVA[index] > std::max(0.52, 0.77-0.025*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 1.479)
    mvaRec = Electron_MVA[index] > std::max(0.11, 0.56-0.045*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 2.5)
    mvaRec = Electron_MVA[index] > std::max(-0.01, 0.48-0.049*(Electron_pt[index]-15));

  return ((Electron_pt[index] > 20) &&
	  (Electron_miniPFRelIso_all[index] < 0.12) &&
	  (mvaRec)  &&
	  (Electron_dz[index] < 0.1) &&
	  (Electron_dxy[index] < 0.05) &&
	  (Electron_sip3d[index] < 4) );
}

bool TTTSelector::IsGoodMVAElectron2017(size_t index) {
  bool mvaRec = false;
  if(abs(Electron_eta[index]) < 0.8)
    mvaRec = std::max(0.52, 0.77-0.025*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 1.479)
    mvaRec = std::max(0.11, 0.56-0.045*(Electron_pt[index]-15));
  else if(abs(Electron_eta[index]) < 2.5)
    mvaRec = std::max(-0.01, 0.48-0.049*(Electron_pt[index]-15));

  return ((Electron_pt[index] > 20) &&
	  mvaRec);
}

//40
bool TTTSelector::IsGoodJet(size_t index) {
  return ((Jet_pt[index] > 40) &&
	  (abs(Jet_eta[index]) < 2.4) &&
	  isOverlap(index)
	  );
}

/// TODO: add toggle for different btag stuff

//25
bool TTTSelector::IsGoodBJet(size_t index) {
  return ((Jet_pt[index] > 25) &&
	  (abs(Jet_eta[index]) < 2.4) &&
	  (Jet_btagCSVV2[index] > 0.8484) &&  
	  //(Jet_btagDeepB[index] > 0.6324) &&
	  isOverlap(index)
	  );
}

bool TTTSelector::isLooseJetId(size_t index) {
  return (Jet_neHEF[index] < 0.99 &&
	  Jet_neEmEF[index] < 0.99 &&
	  Jet_nConstituents[index] > 1 &&
	  Jet_chHEF[index] > 0 &&
	  Jet_chEmEF[index] < 0.99
	  );
}

bool TTTSelector::isTightJetId(size_t index) {
  return ((Jet_neHEF[index] < 0.9) &&
	  (Jet_neEmEF[index] < 0.9) &&
	  (Jet_nConstituents[index] > 1) &&
	  (Jet_chHEF[index] > 0)
	  );
}

bool TTTSelector::passFullIso(LorentzVector& lep, int I2, int I3) {
  LorentzVector closeJet;
  double minDR = 10;
  for(size_t index = 0; index < nJet; index++) {
    LorentzVector jet(Jet_pt[index], Jet_eta[index], Jet_phi[index], Jet_mass[index]);
    double dr = reco::deltaR(jet, lep);
    if(minDR > dr) {
      closeJet = jet;
      minDR = dr;
    }
  }
  
  if(lep.Pt()/closeJet.Pt() > I2 ) return true;
  
  auto diff = closeJet.Vect() - lep.Vect();
  auto cross = diff.Cross(lep.Vect());
  return (cross.Mag2()/diff.Mag2() > I3*I3);
}

bool TTTSelector::isOverlap(size_t index) {
  LorentzVector tmp(Jet_pt[index], Jet_eta[index], Jet_phi[index], Jet_mass[index]);
  double dR = 0.4;
  bool passOverlap = true;
  for(auto lep: goodParts) {
    passOverlap = (reco::deltaR(tmp, lep.v) > dR);
    if(!passOverlap) return false;
  }
  return true;
}

void TTTSelector::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) {
  
  int step = 0;
  Fill1D("CutFlow", 0);

  Fill1D("nGoodParts", goodParts.size());
  
  /// 2 good leptons
  if(goodParts.size() != 2) return;
  Fill1D("CutFlow", ++step);

  
  // std::cout<<"charge array : " <<charge[0]<<""<< charge[1]<<std::endl;
  if(goodParts[0].Charge() * goodParts[1].Charge() >= 0) 
    Fill1D("before", 0);
   
  if(goodParts[0].Charge() * goodParts[1].Charge() <= 0) 
    Fill1D("before", 1);
   
  // first lep requirement 
  if(goodParts[0].Pt() < 25) return;
  Fill1D("CutFlow", ++step);
  
  // same sign requirement
  // if(goodParts[0].Charge() * goodParts[1].Charge() <= 0) return;
  Fill1D("CutFlow", ++step);

  // met cut 50
  if (MET < 50) return;
  Fill1D("CutFlow", ++step);

  // ht cut 
  if(HT <300 ) return;
  Fill1D("CutFlow", ++step);
  
  // jet cut 4
  if(nTightJet < 4) return;
  Fill1D("CutFlow", ++step);
  
  // bjet cut 2
  if(nBJets < 2) return;
  Fill1D("CutFlow", ++step);

  if(goodParts[0].Charge() * goodParts[1].Charge() >= 0)
    Fill1D("after", 0);
  if(goodParts[0].Charge() * goodParts[1].Charge() <=0)
    Fill1D("after", 1);
 
  // // veto cut
  // if(!passesLeptonVeto)
  //   Fill1D("CutFlow", ++step);
  // in SR stuff
  
  Fill1D("ptl1", goodParts[0].Pt());
  Fill1D("ptl2", goodParts[1].Pt());
  Fill1D("SR", getSRBin());
  Fill1D("njet", nTightJet);
  Fill1D("nbjet", nBJets);
  Fill1D("phil1", goodParts[0].Phi());
  Fill1D("phil2", goodParts[1].Phi());
  Fill1D("HT", HT );
  Fill1D("MET", MET);
  Fill1D("nelec", nElectron);
  Fill1D("nmuon", nMuon);
  Fill1D("lept_charge", charge[0]);
  Fill1D("lept_charge", charge[1]);

  Fill2D("bJetvsJets", nTightJet, nBJets);
  

  for(size_t i = 0; i < nJet; i++) {
    if(IsGoodJet(i)) {
      Fill1D("jetpt", Jet_pt[i]);
      Fill1D("jetphi",Jet_phi[i]);
    }
    if(IsGoodBJet(i)) {
      Fill1D("bjetpt", Jet_pt[i]);
      Fill1D("bjetphi",Jet_phi[i]);
    }
  }
 }

void TTTSelector::SetupNewDirectory() {
  SelectorBase::SetupNewDirectory();

  InitializeHistogramsFromConfig();
}

int TTTSelector::getSRBin() {
  if(nBJets == 2) {
    if(nTightJet <= 5)      return 0;
    if(nTightJet == 6)      return 1;
    else if(nTightJet == 7) return 2;
    else if(nTightJet >= 8) return 3;
  } else if(nBJets == 3) {
    if(nTightJet == 5)      return 4;
    else if(nTightJet == 6) return 4;
    else if(nTightJet == 7) return 5;
    else if(nTightJet >= 8) return 5;
  } else if(nBJets >= 4) {
    if(nTightJet >= 5)     return 6;
  }
  return -1;

  // if(nBJets == 2) {
  //   if(nTightJet <= 5)      return 0;
  //   if(nTightJet == 6)      return 1;
  //   else if(nTightJet == 7) return 2;
  //   else if(nTightJet >= 8) return 3;
  // } else if(nBJets == 3) {
  //   if(nTightJet == 5)      return 4;
  //   else if(nTightJet == 6) return 5;
  //   else if(nTightJet == 7) return 6;
  //   else if(nTightJet >= 8) return 7;
  // } else if(nBJets >= 4) {
  //   if(nTightJet >= 5)     return 8;
  // }
  // return -1;

  
 }

