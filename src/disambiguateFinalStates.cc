/*
 * TSelector deduplicate events in combinatorical ntuple
 * by choosing event with l1_l2_Mass pairing closest to
 * M_PDG(Z) = 91.1876
 *
 * Modified from N. Smith, U. Wisconsin
 *
 * https://github.com/nsmith-/ZHinvAnalysis/blob/master/disambiguateFinalStates.C
 */
#include "Analysis/SelectorTools/interface/disambiguateFinalStates.h"

void disambiguateFinalStates::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;

  fChain->SetBranchAddress(zCand_name, &Mass, &b_Mass);
  fChain->SetBranchAddress("evt", &evt, &b_evt);
  fChain->SetBranchAddress("run", &run, &b_run);

  SafeDelete(fCutFormula);
  fCutFormula = new TTreeFormula("CutFormula", fOption, fChain);
  fCutFormula->SetQuickLoad(kTRUE);
  if (!fCutFormula->GetNdim()) { delete fCutFormula; fCutFormula = 0; }
}

Bool_t disambiguateFinalStates::Notify()
{
  return kTRUE;
}

void disambiguateFinalStates::Begin(TTree * /*tree*/)
{
}

void disambiguateFinalStates::SlaveBegin(TTree * /*tree*/)
{
  fBestCandidateEntryList = new TEntryList("bestCandidates", "Entry List of disambiguated combinatoric candidates");
  fOutput->Add(fBestCandidateEntryList);
}

Bool_t disambiguateFinalStates::Process(Long64_t entry)
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
    b_Mass->GetEntry(entry);
    Float_t discriminant = fabs(Mass-91.1876);
    fEntriesToCompare.push_back(entry);
    fEntryDiscriminants.push_back(discriminant);
  }

  if ( entry == fChain->GetEntries()-1 ) {
    findBestEntry();
  }

  return kTRUE;
}

void disambiguateFinalStates::SlaveTerminate()
{
  fBestCandidateEntryList->OptimizeStorage();
  // Pointer is owned by fOutput, dereference
  fBestCandidateEntryList = nullptr;
}

void disambiguateFinalStates::Terminate()
{
}

void disambiguateFinalStates::findBestEntry()
{
  Long64_t bestEntry = -1L;
  Float_t lowestDiscriminant = 1e100;
  for (size_t i=0; i<fEntriesToCompare.size(); ++i)
  {
    if ( lowestDiscriminant > fEntryDiscriminants[i] )
    {
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
}
