/*
 *
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
#ifndef disambiguateFinalStatesZZ_h
#define disambiguateFinalStatesZZ_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TEntryList.h>
#include <TTreeFormula.h>

class disambiguateFinalStatesZZ : public TSelector {
public :
  TTree          *fChain;

  float           l1_l2_Mass;
  float           l1_Pt;
  float           l2_Pt;
  float           l3_l4_Mass;
  float           l3_Pt;
  float           l4_Pt;
  ULong64_t        evt;
  UInt_t           run;

  TBranch        *b_l1_l2_Mass;
  TBranch        *b_l1_Pt; 
  TBranch        *b_l2_Pt;
  TBranch        *b_l3_l4_Mass;
  TBranch        *b_l3_Pt;
  TBranch        *b_l4_Pt;

  TBranch        *b_evt;
  TBranch        *b_run;

  TEntryList     *fBestCandidateEntryList;
  
  // TODO Add TTreeFormula, fix issue when using chain.Process 
  // over multiple files (currently segfaults)
  TTreeFormula   *fCutFormula;


  //disambiguateFinalStatesZZ(TTree * /*tree*/ =0) : fChain(0), fBestCandidateEntryList(0), /*fCutFormula(0),*/ fCurrentRun(-1), fCurrentEvt(-1) { }
  disambiguateFinalStatesZZ(TTree * /*tree*/ =0) : fChain(0), fBestCandidateEntryList(0), fCutFormula(0), fCurrentRun(-1), fCurrentEvt(-1) { }
  //virtual ~disambiguateFinalStatesZZ() { SafeDelete(fBestCandidateEntryList); /*SafeDelete(fCutFormula)*/; }
  virtual ~disambiguateFinalStatesZZ() { SafeDelete(fBestCandidateEntryList); SafeDelete(fCutFormula); }
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
  void setZCandidateBranchName(const char* name1,const char* name2, const char* name3,const char* name4,const char* name5,const char* name6) { 
    l1_l2_Cand_mass = const_cast<char*>(name1);  
    l1_Cand_pt = const_cast<char*>(name2);  
    l2_Cand_pt = const_cast<char*>(name3); 
    l3_l4_Cand_mass = const_cast<char*>(name4); 
    l3_Cand_pt = const_cast<char*>(name5);
    l4_Cand_pt = const_cast<char*>(name6); 
  }

  ClassDef(disambiguateFinalStatesZZ,0);

private :
  void findBestEntry();
  UInt_t fCurrentRun;
  ULong64_t fCurrentEvt;
  std::vector<Long64_t> fEntriesToCompare;
  std::vector<float> fEntryDiscriminants; 
  std::vector<float> fEntryZ2PtSum;
  char* l1_l2_Cand_mass; 
  char* l1_Cand_pt;
  char* l2_Cand_pt; 
  char* l3_l4_Cand_mass; 
  char* l3_Cand_pt;
  char* l4_Cand_pt; 
};

#endif

