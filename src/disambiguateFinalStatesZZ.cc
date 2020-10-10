/*
 * TSelector deduplicate events in combinatorical ntuple
 * by choosing event with l1_l2_Mass pairing closest to
 * M_PDG(Z) = 91.1876 
 * with the highest scalar Pt sum of the remaining leptons
 * used as a tiebreaker. 
 *
 * Modified from K. Long, U. Wisconsin
 *
 * https://https://github.com/kdlong/SelectorTools/blob/master/src/disambiguateFinalStatesZZ.cc
 */
#include "Analysis/SelectorTools/interface/disambiguateFinalStatesZZ.h"

void disambiguateFinalStatesZZ::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;

  fChain->SetBranchAddress(l1_l2_Cand_mass, &l1_l2_Mass, &b_l1_l2_Mass);
  fChain->SetBranchAddress(l1_Cand_pt, &l1_Pt, &b_l1_Pt);
  fChain->SetBranchAddress(l2_Cand_pt, &l2_Pt, &b_l2_Pt);
  fChain->SetBranchAddress(l3_l4_Cand_mass, &l3_l4_Mass, &b_l3_l4_Mass);
  fChain->SetBranchAddress(l3_Cand_pt, &l3_Pt, &b_l3_Pt); 
  fChain->SetBranchAddress(l4_Cand_pt, &l4_Pt, &b_l4_Pt);
  fChain->SetBranchAddress("evt", &evt, &b_evt);
  fChain->SetBranchAddress("run", &run, &b_run);

  SafeDelete(fCutFormula);
  fCutFormula = new TTreeFormula("CutFormula", fOption, fChain);
  fCutFormula->SetQuickLoad(kTRUE);
  if (!fCutFormula->GetNdim()) { delete fCutFormula; fCutFormula = 0; }
}

Bool_t disambiguateFinalStatesZZ::Notify()
{
  return kTRUE;
}

void disambiguateFinalStatesZZ::Begin(TTree * /*tree*/)
{
}

void disambiguateFinalStatesZZ::SlaveBegin(TTree * /*tree*/)
{
  fBestCandidateEntryList = new TEntryList("bestCandidates", "Entry List of disambiguated combinatoric candidates");
  fOutput->Add(fBestCandidateEntryList);
}

Bool_t disambiguateFinalStatesZZ::Process(Long64_t entry)
{
  b_evt->GetEntry(entry);
  b_run->GetEntry(entry);

  if ( !(run == fCurrentRun && evt == fCurrentEvt) )
  {
    findBestEntry();
  }

  fCurrentRun = run;
  fCurrentEvt = evt;
  
  // TODO Understand why this gives segfault for chains
  // with multiple entries
  if ( fCutFormula && fCutFormula->EvalInstance() > 0. )
    {
    b_l1_l2_Mass->GetEntry(entry);
    b_l1_Pt->GetEntry(entry);
    b_l2_Pt->GetEntry(entry);
    b_l3_l4_Mass->GetEntry(entry);
    b_l3_Pt->GetEntry(entry); 
    b_l4_Pt->GetEntry(entry);
    
    float mass_discriminant,Z2ptSum; 
    //This condition identifies the Z1 candidate
    //Required for the 2e2mu state but redundant for the 4e,4mu state however it should be quick comparison
    if(fabs(l1_l2_Mass-91.1876) < fabs(l3_l4_Mass-91.1876)){
      mass_discriminant = fabs(l1_l2_Mass-91.1876);
      Z2ptSum = l3_Pt+l4_Pt;}
    else{ 
      mass_discriminant = fabs(l3_l4_Mass-91.1876);
      Z2ptSum = l1_Pt+l2_Pt;}  
    
    fEntriesToCompare.push_back(entry);
    fEntryDiscriminants.push_back(mass_discriminant);
    fEntryZ2PtSum.push_back(Z2ptSum);
  }

  if ( entry == fChain->GetEntries()-1 ) {
    findBestEntry();
  }

  return kTRUE;
}

void disambiguateFinalStatesZZ::SlaveTerminate()
{
  fBestCandidateEntryList->OptimizeStorage();
  // Pointer is owned by fOutput, dereference
  fBestCandidateEntryList = nullptr;
}

void disambiguateFinalStatesZZ::Terminate()
{
}

void disambiguateFinalStatesZZ::findBestEntry()
{
  //The correct row is the one with Z1 closest
  //to on-shell, with the highest scalar Pt sum of the remaining leptons
  // used as a tiebreaker. 
  Long64_t bestEntry = -1L;
  float lowestDiscriminant = 1e100;
  float MaxPtSum = 0.0;
  for (size_t i=0; i<fEntriesToCompare.size(); ++i)
  {
    if ((fEntryDiscriminants[i] < lowestDiscriminant) || ((fEntryDiscriminants[i] == lowestDiscriminant) && (fEntryZ2PtSum[i] > MaxPtSum)))
    {
      MaxPtSum = fEntryZ2PtSum[i];
      lowestDiscriminant = fEntryDiscriminants[i];
      bestEntry = fEntriesToCompare[i];
    }
  }

  if ( bestEntry >= 0 )
  {
    fBestCandidateEntryList->Enter(bestEntry);
  }

  fEntriesToCompare.clear();
  fEntryDiscriminants.clear();
  fEntryZ2PtSum.clear();
}
