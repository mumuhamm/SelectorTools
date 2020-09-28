
#include "Analysis/SelectorTools/interface/Efficiency.h"
#include "TLorentzVector.h"
#define Fill1D(NAME, VALUE_) HistFullFill(histMap1D_, NAME, variation.first, VALUE_, weight);

#define Fill2D(NAME, VALUE1_, VALUE2_) HistFullFill(histMap2D_, NAME, variation.first, VALUE1_, VALUE2_, weight);

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>LorentzVector;



void Efficiency::Init(TTree *tree) {

    hists1D_ = {"GenElecPt", "GenMuonPt", "MatchElecPt", "MatchMuonPt",
		"FakeElecPt", "FakeMuonPt", "NGenElec", "NGenMuon","GenJetPt", "GenBJetPt",
		"MatchJetPt", "MatchBJetPt", "FakeJetPt", "FakeBJetPt", "NGenBJet", "NGenJet",
		"NMatchBJet", "NFakeBJet"};
    
    hists2D_ = {"BEff_b_btag", "BEff_j_btag", "GenbJetsvsJets"};
    allChannels_ = {{all, "all"}};
    
    b.SetTree(tree);    
    SelectorBase::Init(tree);
    
    TList* slaveClassList = (TList*)GetInputList()->Clone();
    slaveClassList->Add(new TNamed("isSlaveClass", "isSlaveClass"));
    TTTAna.SetInputList(slaveClassList);
    TTTAna.Init(tree);
    TTTAna.SetBranches();
}

void Efficiency::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();

    InitializeHistogramsFromConfig();
}

void Efficiency::SetBranchesNanoAOD() {
    b.CleanUp();
    
    b.SetBranch("nGenPart", nGenPart);
    b.SetBranch("GenPart_pt", GenPart_pt);
    b.SetBranch("GenPart_pdgId", GenPart_pdgId);
    b.SetBranch("GenPart_eta", GenPart_eta);
    b.SetBranch("GenPart_phi", GenPart_phi);
    b.SetBranch("GenPart_mass",GenPart_mass);
    b.SetBranch("GenPart_genPartIdxMother", GenPart_mother);

    b.SetBranch("nGenJet", nGenJet);
    b.SetBranch("GenJet_eta", GenJet_eta);
    b.SetBranch("GenJet_phi", GenJet_phi);
    b.SetBranch("GenJet_pt", GenJet_pt);
    b.SetBranch("GenJet_mass", GenJet_mass);
    b.SetBranch("GenJet_partonFlavour", GenJet_partonFlavour);
}

void Efficiency::clearValues() {
    weight = 1;
    Leptons.clear();
    Jets.clear();
}

void Efficiency::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) {
    clearValues();
    b.SetEntry(entry);

    TTTAna.LoadBranchesNanoAOD(entry,variation);
    
    if (nGenJet > N_KEEP_GEN_JET) {
	std::string message = "Found more Gen jets than max read number. \n Found" ;
	message += std::to_string(nGenJet);
	message += " particles\n  --> Max read number was ";
	message += std::to_string(N_KEEP_GEN_JET);
	message += "\nExiting because this can cause problems. Increase N_KEEP_GEN_JET to avoid this error.\n";
	throw std::domain_error(message);
    }

    for (size_t i=0; i< nGenPart; i++) {
	if (std::abs(GenPart_pdgId[i]) != PID_ELECTRON && std::abs(GenPart_pdgId[i]) != PID_MUON) continue;
	if (GenPart_pt[i] < 20 || abs(GenPart_pdgId[GenPart_mother[i]]) != 24) continue;

	Leptons.push_back(GenPart());
	Leptons.back().SetupGen(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i],
				GenPart_pdgId[i]);
    }

    for (size_t i=0; i < nGenJet; i++) {
	if (GenJet_pt[i] < 25) continue;
	
	Jets.push_back(GenPart());
	Jets.back().SetupGen(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i],
			     GenJet_partonFlavour[i]);
    }

    fillReco(Leptons, TTTAna.goodLeptons);
    fillReco(Jets, TTTAna.goodJets);

    
    weight = TTTAna.weight;
    
    channelName_ = "all";
    channel_ = channelMap_[channelName_];
}

void Efficiency::fillReco(std::vector<GenPart>& genList, const std::vector<GoodPart>& recoList) {
    for(auto reco : recoList) {
	bool isMatched = false;
	for(auto& gen : genList) {
	    if (reco::deltaR(gen.gVector(), reco.v)<0.05) {
		gen.SetupReco(reco);
		isMatched = true;
		break;
	    }
	}
	if(!isMatched) {
	    genList.push_back(GenPart());
	    genList.back().SetupReco(reco);
	}
    }
}

void Efficiency::FillHistograms(Long64_t entry, std::pair<Systematic, std::string> variation) {
    int nGenElec = 0;
    int nGenMuon = 0;
    int nGenBJet = 0;
    int nGenJet = 0;
    int nMatchBJet =0;
    int nFakeBJet = 0;

    for (auto part: Leptons) {
	std::string partName = lepNameMap[part.Id()];
	std::string partPt = partName + "Pt";

	if(part.isMatched()) {
	    Fill1D("Match" + partPt, part.gPt());
	} else if(part.isFaked()) {
	    Fill1D("Fake" + partPt, part.rPt());
	    continue;
	}
	Fill1D("Gen" + partPt, part.gPt());

	if (part.Id() == PID_MUON)  nGenMuon++;
	else 	                    nGenElec++;        
    }

    for(auto jet: Jets) {
	// Jets
	if(jet.isMatched()) {
	    Fill1D("MatchJetPt", jet.gPt());
	} else if(jet.isFaked()) {
	    Fill1D("FakeJetPt", jet.rPt());
	}
	if(!jet.isFaked()) {
	    Fill1D("GenJetPt", jet.gPt());
            nGenJet++;
	}

	// BJet
        if(jet.isBMatched()) {
	    Fill1D("MatchBJetPt", jet.gPt());
	} else if(jet.isBFaked()) {
	    Fill1D("FakeBJetPt", jet.rPt());
	}
        if(jet.isGenB()) {
	    Fill1D("GenBJetPt", jet.gPt());
	    nGenBJet++;
	}
        
    }

    Fill2D("GenbJetsvsJets", nGenJet, nGenBJet);
    Fill1D("NMatchBJet", nMatchBJet);
    Fill1D("NGenBJet", nGenBJet);
    Fill1D("NGenJet", nGenJet);
    Fill1D("NGenElec", nGenElec);
    Fill1D("NGenMuon", nGenMuon);
    Fill1D("NFakeBJet", nFakeBJet);
    return;
}



