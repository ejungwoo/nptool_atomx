void draw_simulation()
{
    bool draw_vertex = true;
    bool draw_kinematics = true;
    bool draw_z_energy = true;
    //bool draw_sep = false;

    //auto file = new TFile("data_simulation/ATOMX_simulation.root");
    auto file = new TFile("data_simulation/simulation_ATOMX.0.root");
    auto tree = (TTree*) file -> Get("SimulatedTree");

    LKDrawingGroup top("draw_simulation");

    TString particle_names[] = {"34Ar", "proton", "37K"};

    if (draw_kinematics)
    {
        for (auto i : {1,2})
        {
            auto hist_tta = new TH1D(Form("hist_tta%d",i),Form("%s;theta",particle_names[i].Data()) ,100, 0, 90);
            if (i==0) tree -> Draw(Form("ReactionConditions.fRC_ThetaCM>>hist_tta%d",i),"","goff");
            else      tree -> Draw(Form("ReactionConditions.fRC_Theta[%d]>>hist_tta%d",i-1,i),"","goff");
            top.AddHist(hist_tta);
        }
        for (auto i : {1,2})
        {
            auto hist_energy = new TH1D(Form("hist_energy%d",i),Form("%s;energy",particle_names[i].Data()), 100, 0, 100);
            if (i==0) tree -> Draw(Form("ReactionConditions.fRC_Beam_Reaction_Energy>>hist_energy%d",i),"","goff");
            else      tree -> Draw(Form("ReactionConditions.fRC_Kinetic_Energy[%d]>>hist_energy%d",i-1,i),"","goff");
            top.AddHist(hist_energy);
        }
        for (auto i : {1,2})
        {
            auto hist_etta = new TH2D(Form("hist_etta%d",i),Form("%s;theta;energy",particle_names[i].Data()), 100, 0, 90, 100, 0, 100);
            if (i==0) tree -> Draw(Form("ReactionConditions.fRC_Beam_Reaction_Energy:ReactionConditions.fRC_ThetaCM>>hist_etta%d",i),"","goff");
            else      tree -> Draw(Form("ReactionConditions.fRC_Kinetic_Energy[%d]:ReactionConditions.fRC_Theta[%d]>>hist_etta%d",i-1,i-1,i),"","goff");
            top.AddHist(hist_etta);
        }
    }

    if (draw_vertex)
    {
        auto hist_vxy = new TH2D("hist_vxy",";vx;vy" ,100, -100, 100 ,100, -100, 100);
        tree -> Draw("ReactionConditions.fRC_Vertex_Position_Y:ReactionConditions.fRC_Vertex_Position_X>>hist_vxy","","goff"); top.AddHist(hist_vxy);

        auto hist_vz = new TH1D("hist_vz",";vz" ,200, -400, 400);
        tree -> Draw("ReactionConditions.fRC_Vertex_Position_Z>>hist_vz","","goff"); top.AddHist(hist_vz);
        hist_vz -> SetMinimum(0);

        //auto hist_3d = new TH3D("hist_3d",";z;x;y" ,200, -400, 400, 200, -400, 400, 200, -400, 400);
        //tree -> Draw("ATOMX.fPosition.x():ATOMX.fPosition.y():ATOMX.fPosition.z()>>hist_3d","","goff"); top.AddHist(hist_3d);
    }

    if (draw_z_energy)
    {
        auto hist_bze = new TH2D("hist_bze",";z;beam energy", 100, -400, 400, 100, 0, 100);
        tree -> Draw("ReactionConditions.fRC_Beam_Reaction_Energy:ReactionConditions.fRC_Vertex_Position_Z>>hist_bze","","goff");
        top.AddHist(hist_bze);
    }

    //auto hist_adiff = new TH1D("hist_adiff",";angle diff", 100, 0, 180);
    //tree -> Draw("(ReactionConditions.GetParticleDirection(0).Angle(ReactionConditions.GetParticleDirection(1))>>hist_adiff","","goff");
    //tree -> Draw("ReactionConditions.GetParticleDirection(0).GetTheta>>hist_adiff","","goff");
    //top.AddHist(hist_adiff);

    top.SetStyle("style");
    top.SetCanvasDivision(4,3);
    //top.SetCanvasDivision(2,5);
    top.Draw();
    top.GetCanvas() -> SaveAs("figures/simulation.png");

    //new TCanvas();
    //tree -> Draw("(ReactionConditions.GetParticleDirection(0)).Angle(ReactionConditions.GetParticleDirection(1))");
    //tree -> Draw("(ReactionConditions.GetParticleDirection(1)).Angle(TVector3(0,0,1))");
    //tree -> Draw("(GetParticleDirection(0)).Angle((GetParticleDirection(1)))");
}
