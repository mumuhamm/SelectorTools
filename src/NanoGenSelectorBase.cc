#include "Analysis/VVAnalysis/interface/NanoGenSelectorBase.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TStyle.h>
#include <regex>

void NanoGenSelectorBase::Init(TTree *tree)
{
    b.SetTree(tree);
    SelectorBase::Init(tree);
    edm::FileInPath mc2hessianCSV("PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv");
    doMC2H_ = name_.find("cp5") == std::string::npos && false;
    std::cout << "INFO: Convert MC to Hessian is " << doMC2H_ << std::endl;
    if (doMC2H_)
        pdfweightshelper_.Init(N_LHEPDF_WEIGHTS_, N_MC2HESSIAN_WEIGHTS_, mc2hessianCSV);
    // NNLOPSLike is just a config name for one MiNNLO sample
    if (name_.find("nnlops") != std::string::npos && name_.find("nnlopslike") == std::string::npos) {
        //std::cout << "INFO: Found NNLOPS sample but not applying weight\n";
        nnlops_ = true;
        std::cout << "INFO: NNLOPS sample will be weighted by NNLO weight\n";
    }
    fReader.SetTree(tree);
}

void NanoGenSelectorBase::SetBranchesNanoAOD() {
}

void NanoGenSelectorBase::LoadBranchesNanoAOD(Long64_t entry, std::pair<Systematic, std::string> variation) { 
    weight = 1;
    fReader.SetLocalEntry(entry);

    channel_ = channelMap_[channelName_];

    bornLeptons.clear();
    bornNeutrinos.clear();
    dressedLeptons.clear();
    jets.clear();

    for (size_t i = 0; i < *nGenDressedLepton; i++) {
        LorentzVector vec;
        if (GenDressedLepton_hasTauAnc.At(i)) {
            continue;
        }
        vec.SetPt(GenDressedLepton_pt.At(i));
        vec.SetEta(GenDressedLepton_eta.At(i));
        vec.SetPhi(GenDressedLepton_phi.At(i));
        vec.SetM(GenDressedLepton_mass.At(i));
        int charge = (GenDressedLepton_pdgId.At(i) < 0) ? 1: -1;
        dressedLeptons.emplace_back(reco::GenParticle(charge, vec, reco::Particle::Point(), GenDressedLepton_pdgId.At(i), 1, true));
    } // No need to sort, they're already pt sorted
    
    if (doBareLeptons_ || doBorn_ || doNeutrinos_ || doPhotons_) {
        bareLeptons.clear();
        fsneutrinos.clear();
        std::vector<unsigned int> idsToKeep = {11, 12, 13, 14};
        if (doPhotons_)
            idsToKeep.push_back(22);

        for (size_t i = 0; i < *nGenPart; i++) {
            bool isHardProcess = (GenPart_statusFlags.At(i) >> 7) & 1;
            if ((doBorn_ && !isHardProcess) || GenPart_status.At(i) != 1)
                continue;
            LorentzVector vec;
            if (std::find(idsToKeep.begin(), idsToKeep.end(), std::abs(GenPart_pdgId.At(i))) != idsToKeep.end()) {
                vec.SetPt(GenPart_pt.At(i));
                vec.SetEta(GenPart_eta.At(i));
                vec.SetPhi(GenPart_phi.At(i));
                vec.SetM(GenPart_mass.At(i));
            }
            if (std::abs(GenPart_pdgId.At(i)) == 11 || std::abs(GenPart_pdgId.At(i)) == 13) {
                int charge = (GenPart_pdgId.At(i) < 0) ? 1: -1;
                auto lep = reco::GenParticle(charge, vec, reco::Particle::Point(), GenPart_pdgId.At(i), GenPart_status.At(i), true);
                if (isHardProcess)
                    bornLeptons.emplace_back(lep);
                if (GenPart_status.At(i) == 1)
                    bareLeptons.emplace_back(lep);
            }
            else if (std::abs(GenPart_pdgId.At(i)) == 12 || std::abs(GenPart_pdgId.At(i)) == 14) {
                auto neutrino = reco::GenParticle(0, vec, reco::Particle::Point(), GenPart_pdgId.At(i), GenPart_status.At(i), true);
                if (isHardProcess)
                    bornNeutrinos.emplace_back(neutrino);
                if (GenPart_status.At(i) == 1)
                    fsneutrinos.emplace_back(neutrino);
            }
            else if (std::abs(GenPart_pdgId.At(i)) == 22) {
                photons.emplace_back(reco::GenParticle(0, vec, reco::Particle::Point(), GenPart_pdgId.At(i), GenPart_status.At(i), true));
            }
        }
        neutrinos = fsneutrinos;
        
        // Sort descending
        auto compareMaxByPt = [](const reco::GenParticle& a, const reco::GenParticle& b) { return a.pt() > b.pt(); };
        std::sort(bareLeptons.begin(), bareLeptons.end(), compareMaxByPt);
        std::sort(bornLeptons.begin(), bornLeptons.end(), compareMaxByPt);

        // Warning! Only really works for the W
        if (bareLeptons.size() > 0 && doPhotons_) {
            auto& lep = bareLeptons.at(0);
            photons.erase(std::remove_if(photons.begin(), photons.end(), 
                    [lep] (const reco::GenParticle& p) { return reco::deltaR(p, lep) > 0.1; }),
                photons.end()
            );
        }
    }

    leptons = dressedLeptons;

    ht = 0;
    for (size_t i = 0; i < *nGenJet; i++) {
        LorentzVector jet;
        jet.SetPt(GenJet_pt.At(i));
        jet.SetEta(GenJet_eta.At(i));
        jet.SetPhi(GenJet_phi.At(i));
        jet.SetM(GenJet_mass.At(i));
        if (jet.pt() > 30 && !helpers::overlapsCollection(jet, leptons, 0.4, nLeptons_)) {
            ht += jet.pt();
            jets.emplace_back(jet);
        }
    } // No need to sort jets, they're already pt sorted

    genMet.SetPt(*MET_fiducialGenPt);
    genMet.SetPhi(*MET_fiducialGenPhi);
    genMet.SetM(0.);
    genMet.SetEta(0.);

    weight = *genWeight;

    if (nnlops_) {
        weight *= LHEScaleWeight.At(9);
    }
    if (doMC2H_)
        buildHessian2MCSet();

    SetComposite();
}

void NanoGenSelectorBase::buildHessian2MCSet() {
    double pdfWeights[N_LHEPDF_WEIGHTS_];
    for (size_t i = 0; i < N_LHEPDF_WEIGHTS_; i++) {
        pdfWeights[i] = LHEPdfWeight[i];
    }
    pdfweightshelper_.DoMC2Hessian(1., const_cast<const double*>(pdfWeights), LHEHessianPdfWeight);
}

void NanoGenSelectorBase::SetupNewDirectory() {
    SelectorBase::SetupNewDirectory();
    AddObject<TH1D>(mcPdfWeights_, "MCweights", "MC pdf weights", 200, 0, 2);
    AddObject<TH1D>(hesPdfWeights_, "Hesweights", "Hessian pdf weights", 200, 0, 2);
    AddObject<TH1D>(scaleWeights_, "scaleweights", "Scale weights", 200, 0, 2);

    InitializeHistogramsFromConfig();
}
