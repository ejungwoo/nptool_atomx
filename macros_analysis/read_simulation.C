#include "simulation_tree.C"

void read_simulation()
{
    auto file = new TFile("data_simulation/ATOMX_simulation.root");
    auto tree = (TTree*) file -> Get("SimulatedTree");
    //if (gSystem->AccessPathName("simulation_tree.h")) {
    //    tree -> MakeClass("simulation_tree");
    //    return;
    //}

    simulation_tree st(tree);
    Long64_t nentries = st.fChain->GetEntriesFast();
    nentries = 4;
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {
        Long64_t ientry = st.LoadTree(jentry);
        if (ientry < 0) break;

        st.fChain -> GetEntry(jentry);
        auto numParticles = st.TrackInfo -> GetParticleMultiplicity();
        cout << ientry << " " << numParticles << endl;
        for (auto iParticle=0; iParticle<numParticles; ++iParticle)
        {
            TString name = st.TrackInfo -> fTI_Particle_Name[iParticle];
            TString producedVolume = st.TrackInfo -> fTI_Volume_Name[iParticle];
            double kin    = st.TrackInfo -> fTI_Kinetic_Energy[iParticle];
            double theta  = st.TrackInfo -> fTI_Theta[iParticle];
            double phi    = st.TrackInfo -> fTI_Phi[iParticle];
            double mass   = st.TrackInfo -> fTI_Mass[iParticle];
            double charge = st.TrackInfo -> fTI_Charge[iParticle];
            double z      = st.TrackInfo -> fTI_Z[iParticle];
            double a      = st.TrackInfo -> fTI_A[iParticle];
            double brho   = st.TrackInfo -> fTI_Brho[iParticle];
            TVector3 mom  = st.TrackInfo -> fTI_Momentum[iParticle];
            TVector3 pos(st.TrackInfo -> fTI_PositionX[iParticle],
                         st.TrackInfo -> fTI_PositionY[iParticle],
                         st.TrackInfo -> fTI_PositionZ[iParticle]);

            //hist -> Fill(pos.x(),pos.y()) << endl;
            cout << "   " << name << " " << producedVolume << " " << kin << " " << brho << endl;
            //mom.Print();
            //pos.Print();
        }

        /*

        vector<int> fTI_Index;

        // emmitted particle
        int    GetParticleMultiplicity() const { return fTI_Kinetic_Energy.size(); }
        string GetParticleName(const int& i) const { return fTI_Particle_Name[i]; }
        double GetKineticEnergy(const int& i) const { return fTI_Kinetic_Energy[i]; }
        TVector3 GetMomentum(const int& i) const { return fTI_Momentum[i]; }
        double   GetMomentumX(const int& i) const { return fTI_Momentum_X[i]; }
        double   GetMomentumY(const int& i) const { return fTI_Momentum_Y[i]; }
        double   GetMomentumZ(const int& i) const { return fTI_Momentum_Z[i]; }

        TVector3 GetParticleDirection(const int& i) const;

        double GetThetaLab_WorldFrame(const int& i) const {
            return (GetParticleDirection(i).Angle(TVector3(0, 0, 1))) / deg;
        }

        unsigned int GetEmittedMult() const { return fTI_Particle_Name.size(); }
        */


    }
}
