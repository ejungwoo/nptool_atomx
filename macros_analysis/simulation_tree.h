//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 17 11:42:25 2025 by ROOT version 6.22/06
// from TTree SimulatedTree/Data created / analysed with the nptool package
// found on file: data_simulation/ATOMX_simulation.root
//////////////////////////////////////////////////////////

#ifndef simulation_tree_h
#define simulation_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TInteractionCoordinates.h"
#include "TReactionConditions.h"
#include "MEventReduced.h"
#include "TObject.h"
#include "MEventReduced.h"
#include "TInitialConditions.h"
#include "TTrackInfo.h"

class simulation_tree {
    public :
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Fixed size dimensions of array or collections stored in the TTree if any.
        static constexpr Int_t kMaxCoboAsad = 3;

        // Declaration of leaf types
        TInteractionCoordinates *InteractionCoordinates;
        TReactionConditions *ReactionConditions;
        //MEventReduced   *data;
        UInt_t          fUniqueID;
        UInt_t          fBits;
        ULong_t         event;
        ULong_t         timestamp;
        Int_t           CoboAsad_;
        UInt_t          CoboAsad_fUniqueID[kMaxCoboAsad];   //[CoboAsad_]
        UInt_t          CoboAsad_fBits[kMaxCoboAsad];   //[CoboAsad_]
        UShort_t        CoboAsad_globalchannelid[kMaxCoboAsad];   //[CoboAsad_]
        Bool_t          CoboAsad_hasSaturation[kMaxCoboAsad];   //[CoboAsad_]
        vector<float>   CoboAsad_peakheight[kMaxCoboAsad];
        vector<float>   CoboAsad_peaktime[kMaxCoboAsad];
        TInitialConditions *InitialConditions;
        TTrackInfo      *TrackInfo;
        Int_t           Run;

        // List of branches
        TBranch        *b_InteractionCoordinates;   //!
        TBranch        *b_ReactionConditions;   //!
        TBranch        *b_data_fUniqueID;   //!
        TBranch        *b_data_fBits;   //!
        TBranch        *b_data_event;   //!
        TBranch        *b_data_timestamp;   //!
        TBranch        *b_data_CoboAsad_;   //!
        TBranch        *b_CoboAsad_fUniqueID;   //!
        TBranch        *b_CoboAsad_fBits;   //!
        TBranch        *b_CoboAsad_globalchannelid;   //!
        TBranch        *b_CoboAsad_hasSaturation;   //!
        TBranch        *b_CoboAsad_peakheight;   //!
        TBranch        *b_CoboAsad_peaktime;   //!
        TBranch        *b_InitialConditions;   //!
        TBranch        *b_TrackInfo;   //!
        TBranch        *b_Run;   //!

        simulation_tree(TTree *tree=0);
        virtual ~simulation_tree();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef simulation_tree_cxx
simulation_tree::simulation_tree(TTree *tree) : fChain(0) 
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_simulation/ATOMX_simulation.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("data_simulation/ATOMX_simulation.root");
        }
        f->GetObject("SimulatedTree",tree);

    }
    Init(tree);
}

simulation_tree::~simulation_tree()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t simulation_tree::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t simulation_tree::LoadTree(Long64_t entry)
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

void simulation_tree::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    InteractionCoordinates = 0;
    ReactionConditions = 0;
    InitialConditions = 0;
    TrackInfo = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("InteractionCoordinates", &InteractionCoordinates, &b_InteractionCoordinates);
    fChain->SetBranchAddress("ReactionConditions", &ReactionConditions, &b_ReactionConditions);
    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_data_fUniqueID);
    fChain->SetBranchAddress("fBits", &fBits, &b_data_fBits);
    fChain->SetBranchAddress("event", &event, &b_data_event);
    fChain->SetBranchAddress("timestamp", &timestamp, &b_data_timestamp);
    fChain->SetBranchAddress("CoboAsad", &CoboAsad_, &b_data_CoboAsad_);
    fChain->SetBranchAddress("CoboAsad.fUniqueID", CoboAsad_fUniqueID, &b_CoboAsad_fUniqueID);
    fChain->SetBranchAddress("CoboAsad.fBits", CoboAsad_fBits, &b_CoboAsad_fBits);
    fChain->SetBranchAddress("CoboAsad.globalchannelid", CoboAsad_globalchannelid, &b_CoboAsad_globalchannelid);
    fChain->SetBranchAddress("CoboAsad.hasSaturation", CoboAsad_hasSaturation, &b_CoboAsad_hasSaturation);
    fChain->SetBranchAddress("CoboAsad.peakheight", CoboAsad_peakheight, &b_CoboAsad_peakheight);
    fChain->SetBranchAddress("CoboAsad.peaktime", CoboAsad_peaktime, &b_CoboAsad_peaktime);
    fChain->SetBranchAddress("InitialConditions", &InitialConditions, &b_InitialConditions);
    fChain->SetBranchAddress("TrackInfo", &TrackInfo, &b_TrackInfo);
    fChain->SetBranchAddress("Run", &Run, &b_Run);
    Notify();
}

Bool_t simulation_tree::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void simulation_tree::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t simulation_tree::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef simulation_tree_cxx
