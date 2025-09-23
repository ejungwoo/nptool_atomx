#include "ASimulatedTree.C"

void read_simulation(TString inputFileName="data/ATOMX_34Ar.sim.root")
{
    bool test = false;

    auto file = new TFile(inputFileName);
    auto tree = (TTree*) file -> Get("SimulatedTree");
    auto st = new ASimulatedTree(tree);
    auto numEvents = tree -> GetEntries();

    if (test) {
        cout << "running read_simulation.C(\""<< inputFileName << "\") in test mode" << endl;
        numEvents = 1;
    }

    auto histZX = new TH2D("hist_atomx_zx",";z (mm);x (mm)",200,-250,250,200,-250,250);
    auto histZY = new TH2D("hist_atomx_zy",";z (mm);y (mm)",200,-250,250,200,-200,200);
    auto histELoss = new TH1D("hist_atomx_eloss",";energy loss (MeV)",200,0,0.2);
    auto histZELoss = new TH2D("hist_atomx_zeloss",";z (mm);energy loss (MeV)",200,-250,250,200,0,0.2);
    auto FillATOMX = [histZY,histZX,histELoss,histZELoss](double energy, double time, TVector3 position) {
        histZX -> Fill(position.z(),position.x());
        histZY -> Fill(position.z(),position.y());
        if (energy>0) {
            histELoss -> Fill(energy);
            histZELoss -> Fill(position.z(),energy);
        }
    };

    for (auto iEvent=0; iEvent<numEvents; ++iEvent)
    {
        st -> GetEntry(iEvent);

        auto numTracks = st -> TrackInfo -> GetParticleMultiplicity();
        auto numSTARK  = st -> STARK -> GetMult();
        auto numATOMX = st -> ATOMX -> GetEnergyLossMultiplicity();
        if (test) {
            cout << "Event " << iEvent << endl;
            cout << "numTracks = " << numTracks << endl;
            cout << "numSTARK  = " << numSTARK  << endl;
            cout << "numATOMX  = " << numATOMX  << endl;
            cout << endl; st -> ReactionConditions -> Dump();
            cout << endl; st -> TrackInfo -> Dump();
            cout << endl; st -> STARK -> Dump();
            cout << endl; st -> ATOMX -> Dump(20);
        }

        for (auto iTrack=0; iTrack<numTracks; ++iTrack)
        {
            TString name = st->TrackInfo -> fTI_Particle_Name[iTrack];
            auto ke      = st->TrackInfo -> fTI_Kinetic_Energy[iTrack];
            auto a       = st->TrackInfo -> fTI_A[iTrack];
            auto z       = st->TrackInfo -> fTI_Z[iTrack];
            auto tta     = st->TrackInfo -> fTI_Theta[iTrack];
            auto phi     = st->TrackInfo -> fTI_Phi[iTrack];
            //TODO
        }
        for (auto iATOMX=0; iATOMX<numATOMX; ++iATOMX)
        {
            double energy     = st->ATOMX -> fEnergyLoss[iATOMX];
            double time       = st->ATOMX -> fTime[iATOMX];
            TVector3 position = st->ATOMX -> fPosition[iATOMX];
            FillATOMX(energy,time,position);
            //TODO
        }
        for (auto iSTARK=0; iSTARK<numSTARK; ++iSTARK)
        {
            int    type = st->STARK -> GetType(iSTARK);
            int    detn = st->STARK -> GetDetN(iSTARK);
            int    fstn = st->STARK -> GetFStN(iSTARK);
            int    bstn = st->STARK -> GetBStN(iSTARK);
            double fre  = st->STARK -> GetFrE (iSTARK);
            double bke  = st->STARK -> GetBkE (iSTARK);
            double upe  = st->STARK -> GetUpE (iSTARK);
            double dwe  = st->STARK -> GetDwE (iSTARK);
            TVector3 pos= st->STARK -> GetPos (iSTARK);
            //TODO
        }
    }

    auto cvs = new TCanvas("cvs","",1500,1000);
    cvs -> Divide(2,2);
    cvs -> cd(1)->SetLogz(); histZX -> Draw("colz");
    cvs -> cd(2)->SetLogz(); histZY -> Draw("colz");
    cvs -> cd(3); histELoss -> Draw();
    cvs -> cd(4)->SetLogz(); histZELoss -> Draw("colz");
    cvs -> SaveAs("figures/ATOMX_34Ar_read_simulation_example.png");
}
