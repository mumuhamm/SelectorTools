#include "Analysis/VVAnalysis/interface/NanoGenSelectorBase.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TParameter.h"
#include <TStyle.h>
#include <regex>

void NanoGenSelectorBase::Init(TTree *tree)
{
    refWeight = 1;
    TParameter<bool>* doTheory = (TParameter<bool>*) GetInputList()->FindObject("theoryUnc");
    doTheoryVars_ = doTheory != nullptr && doTheory->GetVal();
    if (doTheoryVars_) {
        theoryVarSysts_.insert(theoryVarSysts_.begin(), Central);
        theoryVarSysts_.insert(theoryVarSysts_.end(), LHEParticles);
        theoryVarSysts_.insert(theoryVarSysts_.end(), BareLeptons);
    }

    TParameter<bool>* doPtVSplit = (TParameter<bool>*) GetInputList()->FindObject("theoryPtV");
    if (doTheoryVars_ && doPtVSplit != nullptr && doPtVSplit->GetVal()) {
        // For now it's based on the LHE kinematics
        TNamed* selection = (TNamed *) GetInputList()->FindObject("selection");
        std::string selectionName = selection == nullptr ? "" : selection->GetTitle();
        SystMap ptvars = {
            {ptV0to3, "ptV0to3"},
            {ptV3to5, "ptV3to5"},
            {ptV5to7, "ptV5to7"},
            {ptV7to9, "ptV7to9"},
            {ptV9to12, "ptV9to12"},
            {ptV12to15, "ptV12to15"},
            {ptV15to20, "ptV15to20"},
            {ptV20to27, "ptV20to27"},
            {ptV27to40, "ptV27to40"},
            {ptV40toInf, "ptV40toInf"},
        };
        if (selectionName == "None") {
            ptvars = {
                {ptV0to3_lhe, "lhe_ptV0to3"},
                {ptV3to5_lhe, "lhe_ptV3to5"},
                {ptV5to7_lhe, "lhe_ptV5to7"},
                {ptV7to9_lhe, "lhe_ptV7to9"},
                {ptV9to12_lhe, "lhe_ptV9to12"},
                {ptV12to15_lhe, "lhe_ptV12to15"},
                {ptV15to20_lhe, "lhe_ptV15to20"},
                {ptV20to27_lhe, "lhe_ptV20to27"},
                {ptV27to40_lhe, "lhe_ptV27to40"},
                {ptV40toInf_lhe, "lhe_ptV40toInf"},
            };
        }
        for (auto& var : ptvars) {
            systematics_[var.first] = var.second;
            theoryVarSysts_.insert(theoryVarSysts_.end(), var.first);
        }
    }
    b.SetTree(tree);
    altScaleWeights_ = (tree->GetListOfBranches()->FindObject("nLHEScaleWeightAltSet1") != nullptr);
    unknownWeights_ = (tree->GetListOfBranches()->FindObject("nLHEUnknownWeight") != nullptr);
    unknownWeightsAlt_ = (tree->GetListOfBranches()->FindObject("nLHEUnknownWeightAltSet1") != nullptr);
    pdfWeights_ = (tree->GetListOfBranches()->FindObject("nLHEPdfWeight") != nullptr);

    SelectorBase::Init(tree);

    edm::FileInPath mc2hessianCSV("PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv");
    // Off for now
    doMC2H_ = name_.find("cp5") == std::string::npos && false;
    if (doMC2H_) {
        std::cout << "INFO: Will convert MC PDF set to Hessian with MC2Hessian\n";
        pdfweightshelper_.Init(N_LHEPDF_WEIGHTS_, N_MC2HESSIAN_WEIGHTS_, mc2hessianCSV);
    }
    centralWeightIndex_ = -1;
    // NNLOPSLike is just a config name for one MiNNLO sample
    if (name_.find("nnlops") != std::string::npos && name_.find("nnlopslike") == std::string::npos) {
        if (name_.find("nloOnly") == std::string::npos) {
            centralWeightIndex_ = 9;
            std::cout << "INFO: NNLOPS sample will be weighted by NNLO weight\n";
        }
        else
            std::cout << "INFO: Found NNLOPS sample but not applying weight\n";
    }
    else if (name_.find("ref_nnpdf31") != std::string::npos) {
        centralWeightIndex_ = 0;
        std::cout << "INFO: Sample will be weighted by 1st LHE weight\n";
    }

    doFiducial_ = selection_ != None;
    if (!doFiducial_)
        std::cout << "INFO: No fiducial selection will be applied\n";
    //doBorn_ = std::find_if(systematics_.begin(), systematics_.end(), [](auto& s) { return s.first == BornParticles; }) != systematics_.end();
    //doBareLeptons_ = std::find_if(systematics_.begin(), systematics_.end(), [](auto& s) { return s.first == BareLeptons; }) != systematics_.end();
    doBorn_ = true;
    doBareLeptons_ = true;
    fReader.SetTree(tree);
}

void NanoGenSelectorBase::SetBranchesNanoAOD() {
    if (altScaleWeights_) {
        b.SetSpecificBranch("nLHEScaleWeightAltSet1", nLHEScaleWeightAltSet1);
        b.SetSpecificBranch("LHEScaleWeightAltSet1", LHEScaleWeightAltSet1);
    }
    if (pdfWeights_) {
        b.SetSpecificBranch("nLHEPdfWeight", nLHEPdfWeight);
        b.SetSpecificBranch("LHEPdfWeight", LHEPdfWeight);
    }
    if (unknownWeights_) {
        b.SetSpecificBranch("nLHEUnknownWeight", nLHEUnknownWeight);
        b.SetSpecificBranch("LHEUnknownWeight", LHEUnknownWeight);
        if (unknownWeightsAlt_) {
            b.SetSpecificBranch("nLHEUnknownWeightAltSet1", nLHEUnknownWeightAltSet1);
            b.SetSpecificBranch("LHEUnknownWeightAltSet1", LHEUnknownWeightAltSet1);
        }
    }
}

void NanoGenSelectorBase::LoadBranchesNanoAOD(Long64_t entry, SystPair variation) { 
    weight = 1;
    if (altScaleWeights_) {
        b.SetSpecificEntry(entry, "nLHEScaleWeightAltSet1");
        b.SetSpecificEntry(entry, "LHEScaleWeightAltSet1");
    }
    if (pdfWeights_) {
        b.SetSpecificEntry(entry, "LHEPdfWeight");
        b.SetSpecificEntry(entry, "nLHEPdfWeight");
    }
    if (unknownWeights_) {
        b.SetSpecificEntry(entry, "LHEUnknownWeight");
        b.SetSpecificEntry(entry, "nLHEUnknownWeight");
        if (unknownWeightsAlt_) {
            b.SetSpecificEntry(entry, "LHEUnknownWeightAltSet1");
            b.SetSpecificEntry(entry, "nLHEUnknownWeightAltSet1");
        }
    }

    if (!altScaleWeights_)
        nLHEScaleWeightAltSet1 = 0;
    if (!pdfWeights_)
        nLHEPdfWeight = 0;
    if (!unknownWeights_)
        nLHEUnknownWeight = 0;
    if (!unknownWeightsAlt_)
        nLHEUnknownWeightAltSet1 = 0;
    fReader.SetLocalEntry(entry);

    channel_ = channelMap_[channelName_];

    // Sort descending
    auto compareMaxByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() > b.pt(); };
    std::vector<unsigned int> idsToKeep = {11, 12, 13, 14};
    if (doPhotons_)
        idsToKeep.push_back(22);

    // This is a bit nasty because it assumes that Central is the first variation, which is NOT guaranteed!
    if (variation.first == Central) {
        bornLeptons.clear();
        dressedLeptons.clear();
        bareLeptons.clear();

        bornNeutrinos.clear();
        fsneutrinos.clear();

        jets.clear();
        for (size_t i = 0; i < *nGenDressedLepton; i++) {
            LorentzVector vec;
            if (GenDressedLepton_hasTauAnc.At(i)) {
                continue;
            }
            dressedLeptons.emplace_back(makeGenParticle(GenDressedLepton_pdgId.At(i), 1, GenDressedLepton_pt.At(i), 
                                        GenDressedLepton_eta.At(i), GenDressedLepton_phi.At(i), GenDressedLepton_mass.At(i)));
        } // No need to sort, they're already pt sorted
        leptons = dressedLeptons;
        // Do once, only if bare or born leptons are requested
        if (doNeutrinos_ || doBorn_ || doBareLeptons_) {
            for (size_t i = 0; i < *nGenPart; i++) {
                bool isHardProcess = (GenPart_statusFlags.At(i) >> 7) & 1;
                bool isPrompt = (GenPart_statusFlags.At(i) >> 0) & 1;
                if (!(isHardProcess || (GenPart_status.At(i) == 1 && isPrompt)))
                    continue;
                if (std::find(idsToKeep.begin(), idsToKeep.end(), std::abs(GenPart_pdgId.At(i))) == idsToKeep.end())
                    continue;

                auto part = makeGenParticle(GenPart_pdgId.At(i), GenPart_status.At(i), GenPart_pt.At(i), 
                        GenPart_eta.At(i), GenPart_phi.At(i), GenPart_mass.At(i));
                if (std::abs(part.pdgId()) == 11 || std::abs(part.pdgId()) == 13) {
                    if (doBareLeptons_ && isPrompt && GenPart_status.At(i) == 1)
                        bareLeptons.emplace_back(part);
                    if (isHardProcess && doBorn_)
                        bornLeptons.emplace_back(part);
                }
                else if (std::abs(part.pdgId()) == 12 || std::abs(part.pdgId()) == 14) {
                    if (GenPart_status.At(i) == 1)
                        fsneutrinos.emplace_back(part);
                    if (isHardProcess && doBorn_)
                        bornNeutrinos.emplace_back(part);
                }
                else if (std::abs(part.pdgId()) == 22 && isPrompt) {
                    photons.emplace_back(part);
                }
            }
        }
        // Warning! Only really works for the W
        if (bareLeptons.size() > 0 && doPhotons_) {
            auto& lep = bareLeptons.at(0);
            photons.erase(std::remove_if(photons.begin(), photons.end(), 
                    [lep] (const reco::GenParticle& p) { return reco::deltaR(p, lep) > 0.1; }),
                photons.end()
            );
        }
        neutrinos = fsneutrinos;
        std::sort(bareLeptons.begin(), bareLeptons.end(), compareMaxByPt);
        std::sort(bornLeptons.begin(), bornLeptons.end(), compareMaxByPt);
    }
    else if (variation.first == BareLeptons) {
        leptons = bareLeptons;
        neutrinos = fsneutrinos;
    }
    else if (variation.first == BornParticles) {
        leptons = bornLeptons;
        neutrinos = bornNeutrinos;
    }
    else if (variation.first == LHEParticles) {
        lheLeptons.clear();
        lheNeutrinos.clear();
        for (size_t i = 0; i < *nLHEPart; i++) {
            if (std::find(idsToKeep.begin(), idsToKeep.end(), std::abs(LHEPart_pdgId.At(i))) == idsToKeep.end())
                continue;

            auto part = makeGenParticle(LHEPart_pdgId.At(i), 1, LHEPart_pt.At(i), 
                    LHEPart_eta.At(i), LHEPart_phi.At(i), LHEPart_mass.At(i));

            if (std::abs(part.pdgId()) == 11 || std::abs(part.pdgId()) == 13) {
                lheLeptons.emplace_back(part);
            }
            else if (std::abs(part.pdgId()) == 12 || std::abs(part.pdgId()) == 14) {
                lheNeutrinos.emplace_back(part);
            }
            // No dR check here, should add
            else if ((std::abs(part.pdgId()) < 6 || std::abs(part.pdgId()) == 21) && part.pt() > 30) {
                jets.emplace_back(part.polarP4());
            }
        }
        std::sort(lheLeptons.begin(), lheLeptons.end(), compareMaxByPt);
        leptons = lheLeptons;
        neutrinos = lheNeutrinos;
    }
    else if (variation.first == ptV0to3_lhe ||
            variation.first == ptV3to5_lhe || 
            variation.first == ptV5to7_lhe ||
            variation.first == ptV7to9_lhe ||
            variation.first == ptV9to12_lhe ||
            variation.first == ptV12to15_lhe || 
            variation.first == ptV15to20_lhe ||
            variation.first == ptV20to27_lhe ||
            variation.first == ptV27to40_lhe ||
            variation.first == ptV40toInf_lhe ) {
        leptons = lheLeptons;
        neutrinos = lheNeutrinos;
    }
    else {
        leptons = dressedLeptons;
        neutrinos = fsneutrinos;
    }
        
    if (variation.first != LHEParticles) {
        jets.clear();
        ht = 0;
        for (size_t i = 0; i < *nGenJet; i++) {
            auto jet = makeGenParticle(0, 1, GenJet_pt.At(i), 
                    GenJet_eta.At(i), GenJet_phi.At(i), GenJet_mass.At(i));
            if (jet.pt() > 30 && !helpers::overlapsCollection(jet.polarP4(), leptons, 0.4, nLeptons_)) {
                ht += jet.pt();
                jets.emplace_back(jet.polarP4());
            }
        } // No need to sort jets, they're already pt sorted

        genMet.SetPt(*MET_fiducialGenPt);
        genMet.SetPhi(*MET_fiducialGenPhi);
        genMet.SetM(0.);
        genMet.SetEta(0.);
    }

    weight = *genWeight;       
    if (refWeight == 1)
        refWeight = weight;

    if (centralWeightIndex_ != -1) {
        weight *= LHEScaleWeight.At(centralWeightIndex_);
    }
    if (doMC2H_)
        buildHessian2MCSet();

    SetComposite();
}

reco::GenParticle NanoGenSelectorBase::makeGenParticle(int pdgid, int status, float pt, float eta, float phi, float m) {
    LorentzVector vec;
    vec.SetPt(pt);
    vec.SetEta(eta);
    vec.SetPhi(phi);
    vec.SetM(m);
    int charge = (pdgid < 0) ? 1: -1;
    if (std::abs(pdgid) ==  12 || std::abs(pdgid) ==  14 || std::abs(pdgid) ==  16 || std::abs(pdgid) ==  22 || std::abs(pdgid) ==  23)
        charge = 0;
    
    auto lep = reco::GenParticle(charge, vec, reco::Particle::Point(), pdgid, status, true);
    return lep;
}

void NanoGenSelectorBase::buildHessian2MCSet() {
    double pdfWeights[N_LHEPDF_WEIGHTS_];
//    for (size_t i = 0; i < N_LHEPDF_WEIGHTS_; i++) {
//        pdfWeights[i] = LHEPdfWeight[i];
//    }
    pdfweightshelper_.DoMC2Hessian(1., const_cast<const double*>(pdfWeights), LHEHessianPdfWeight);
}

double NanoGenSelectorBase::breitWignerWeight(double offset) {

    double targetMass = MV_GEN_ + offset;
    double s_hat = mVlhe*mVlhe;
    double offshell = s_hat - MV_GEN_*MV_GEN_;
    double offshellOffset = s_hat - targetMass*targetMass;
    double weight = (offshell*offshell + GAMMAV_GEN_*GAMMAV_GEN_*MV_GEN_*MV_GEN_)/
            (offshellOffset*offshellOffset + GAMMAV_GEN_*GAMMAV_GEN_*targetMass*targetMass);
    return weight;
}

void NanoGenSelectorBase::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();
    AddObject<TH1D>(mcWeights_, "genWeights", "gen weights", 200, -10, 10);
    AddObject<TH1D>(mcPdfWeights_, "MCPdfweights", "MC pdf weights", 200, 0, 2);
    AddObject<TH1D>(hesPdfWeights_, "Hesweights", "Hessian pdf weights", 200, 0, 2);
    AddObject<TH1D>(scaleWeights_, "scaleweights", "Scale weights", 200, 0, 2);

    InitializeHistogramsFromConfig();
}
