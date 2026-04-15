//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 17 20:35:03 2026 by ROOT version 6.38.00
// from TTree SimulatedTree/Data created / analysed with the nptool package
// found on file: data/ATOMX_12C.sim.root
//////////////////////////////////////////////////////////

#ifndef SimulatedTree_h
#define SimulatedTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TInteractionCoordinates.h"
#include "TReactionConditions.h"
#include "TATOMXData.h"
#include "TSTARKData.h"
#include "TSTARKRaw.h"
#include "TInitialConditions.h"
#include "TTrackInfo.h"

class SimulatedTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TInteractionCoordinates *InteractionCoordinates;
   TReactionConditions *ReactionConditions;
   TATOMXData      *ATOMX;
   TSTARKData      *STARK;
   TSTARKRaw       *Raw;
   TInitialConditions *InitialConditions;
   Int_t           Run;
   TTrackInfo      *TrackInfo;

   // List of branches
   TBranch        *b_InteractionCoordinates;   //!
   TBranch        *b_ReactionConditions;   //!
   TBranch        *b_ATOMX;   //!
   TBranch        *b_STARK;   //!
   TBranch        *b_Raw;   //!
   TBranch        *b_InitialConditions;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_TrackInfo;   //!

   SimulatedTree(TTree *tree=0);
   virtual ~SimulatedTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SimulatedTree_cxx
SimulatedTree::SimulatedTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data/ATOMX_12C.sim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/ATOMX_12C.sim.root");
      }
      f->GetObject("SimulatedTree",tree);

   }
   Init(tree);
}

SimulatedTree::~SimulatedTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SimulatedTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SimulatedTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SimulatedTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.

   // Set object pointer
   InteractionCoordinates = 0;
   ReactionConditions = 0;
   ATOMX = 0;
   STARK = 0;
   Raw = 0;
   InitialConditions = 0;
   TrackInfo = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("InteractionCoordinates", &InteractionCoordinates, &b_InteractionCoordinates);
   fChain->SetBranchAddress("ReactionConditions", &ReactionConditions, &b_ReactionConditions);
   fChain->SetBranchAddress("ATOMX", &ATOMX, &b_ATOMX);
   fChain->SetBranchAddress("STARK", &STARK, &b_STARK);
   fChain->SetBranchAddress("Raw", &Raw, &b_Raw);
   fChain->SetBranchAddress("InitialConditions", &InitialConditions, &b_InitialConditions);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("TrackInfo", &TrackInfo, &b_TrackInfo);
   Notify();
}

bool SimulatedTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be for a new TTree in a TChain. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void SimulatedTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SimulatedTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SimulatedTree_cxx
