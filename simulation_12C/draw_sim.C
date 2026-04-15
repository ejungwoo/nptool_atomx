#include "SimulatedTree.C"

void draw_sim()
{
    bool verbose = true;
    auto file = new TFile("data/ATOMX_12C.sim.root","read");
    auto tree = (TTree*) file -> Get("SimulatedTree");
    auto t = new SimulatedTree(tree);
    auto numEvents = tree -> GetEntries();
    //if (verbose && numEvents>1) numEvents = 10;

    TATOMXData          *ATOMX;
    TSTARKData          *STARK;
    TSTARKRaw           *Raw;
    TInitialConditions *InitialConditions;
    Int_t               Run;
    TTrackInfo          *TrackInfo;

    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        t->GetEntry(iEvent);
        if (verbose) cout << endl;
        if (verbose) cout << "Event " << iEvent << endl;
        auto b_ti_mult= t->TrackInfo -> GetParticleMultiplicity();
        auto b_ti_ke  = t->TrackInfo -> fTI_Kinetic_Energy[0];
        auto b_ti_a   = t->TrackInfo -> fTI_A[0];
        auto b_ti_z   = t->TrackInfo -> fTI_Z[0];
        TString ti_name  = t->TrackInfo -> fTI_Particle_Name[0];
        auto b_ti_tta = t->TrackInfo -> fTI_Theta[0];
        auto b_ti_phi = t->TrackInfo -> fTI_Phi[0];
        if (verbose) cout << "name=" << ti_name << " a=" << b_ti_a << " z=" << b_ti_z << endl;
        if (verbose) cout << "ke=" << b_ti_ke << " tta=" << b_ti_tta << endl;

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        auto b_multSTARK = t->STARK -> GetMult();
        if (verbose && b_multSTARK==0) continue;
        if (verbose) cout << "mult_STARK=" << b_multSTARK << endl;
        if (verbose) t->STARK -> Dump();
        for (auto iSTARK=0; iSTARK<b_multSTARK; ++iSTARK)
        {
            //int    type = t->STARK -> GetType(iSTARK);
            int    detn = t->STARK -> GetDetN(iSTARK);
            int    fstn = t->STARK -> GetFStN(iSTARK);
            int    bstn = t->STARK -> GetBStN(iSTARK);
            double fre  = t->STARK -> GetFrE (iSTARK);
            double bke  = t->STARK -> GetBkE (iSTARK);
            double upe  = t->STARK -> GetUpE (iSTARK);
            double dwe  = t->STARK -> GetDwE (iSTARK);
            //double t    = t->STARK -> GetT   (iSTARK);
            TVector3 pos = t->STARK -> GetPos (iSTARK);
        }

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        auto b_multATOMX = t->ATOMX -> GetEnergyLossMultiplicity();
        if (verbose) cout << "mult_ATOMX=" << b_multATOMX << endl;
        if (verbose) t->ATOMX -> Dump();
        if (0)
            for (auto iATOMX=0; iATOMX<b_multATOMX; ++iATOMX)
            {
                TVector3 pos    = t->ATOMX -> fPosition[iATOMX];
                double   energy = t->ATOMX -> fEnergyLoss[iATOMX];
                double   time   = t->ATOMX -> fTime[iATOMX];
            }

        if (verbose) break;
    }
}
