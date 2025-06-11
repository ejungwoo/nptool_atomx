#include "simulation_tree.C"

void read_simulation()
{
    bool printReactionConditions = true;
    bool printTrackInfo = true;
    bool printATOMX = true;

    //auto file = new TFile("data_simulation/ATOMX_viewer.root");
    auto file = new TFile(" data_simulation/viewer_ATOMX.0.root");
    auto tree = (TTree*) file -> Get("SimulatedTree");
    //if (gSystem->AccessPathName("simulation_tree.h")) {
    //    tree -> MakeClass("simulation_tree");
    //    return;
    //}

    simulation_tree st(tree);
    auto numEntries = st.fChain->GetEntriesFast();
    numEntries = 10;
    for (Long64_t jEntry=0; jEntry<numEntries; jEntry++)
    {
        Long64_t iEntry = st.LoadTree(jEntry);
        if (iEntry < 0) break;
        st.fChain -> GetEntry(jEntry);
        cout << endl;
        cout << "== Event-" << iEntry << "(" << jEntry << ")" << endl;

        if (printReactionConditions)
        {
            TString name = st.ReactionConditions -> GetBeamParticleName();
            double energy = st.ReactionConditions -> GetBeamEnergy();
            double thetaCM = st.ReactionConditions -> GetThetaCM();
            TVector3 vertex = st.ReactionConditions -> GetVertexPosition();
            cout << "[ReactionConditions]" << endl;
            cout << "  beam name=" << name << ", energy=" << energy << ", thetaCM=" << thetaCM << endl;
            cout << "  vertex=(" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ")" << endl;

            auto numReactionParticles = st.ReactionConditions -> GetParticleMultiplicity();
            for (auto iParticle=0; iParticle<numReactionParticles; ++iParticle)
            {
                name = st.ReactionConditions -> GetParticleName(iParticle);
                double theta = st.ReactionConditions -> GetTheta(iParticle);
                double phi = st.ReactionConditions -> GetPhi(iParticle);
                energy = st.ReactionConditions -> GetKineticEnergy(iParticle);
                cout << "    " << iParticle << ") name=" << name << ", e=" << energy << ", theta=" << theta << ", phi=" << phi << endl;
            }
        }

        if (printTrackInfo)
        {
            auto numParticles = st.TrackInfo -> GetParticleMultiplicity();
            cout << "[TrackInfo]" << endl;
            for (auto iParticle=0; iParticle<numParticles; ++iParticle)
            {
                int id = st.TrackInfo -> fTI_Index[iParticle];
                id = id - iEntry;
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
                cout << "    " << iParticle << ") id=" << id << ", name=" << name << ", volume=" << producedVolume << " e=" << kin << " z=" << pos.z() << endl;
            }
        }

        if (printATOMX)
        {
            auto numPoints = st.ATOMX->GetEnergyLossMultiplicity();
            cout << "[ATOMX]" << endl;
            cout << "  ATOMX containing " << numPoints << " points" << endl;
            auto numReducedPoints = numPoints;
            if (numPoints>20)
                numReducedPoints = 40;
            for (auto iPoint=0; iPoint<numReducedPoints; ++iPoint)
            {
                double eloss = st.ATOMX->fEnergyLoss[iPoint];
                double time = st.ATOMX->fTime[iPoint];
                TVector3 position = st.ATOMX->fPosition[iPoint];
                cout << "    " << iPoint << ") e=" << eloss << ", t=" << time << ", z=" << position.z() << endl;
            }
            if (numReducedPoints!=numPoints)
                cout << "    ..." << endl;
        }
    }
}
