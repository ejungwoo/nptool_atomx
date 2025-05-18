void draw_simulation()
{
    gStyle -> SetPalette(kRainBow);
    gStyle -> SetOptStat(0);

    auto file = new TFile("data_simulation/ATOMX_simulation.root");
    auto tree = (TTree*) file -> Get("SimulatedTree");

    TString select_not_beam = "TrackInfo.fTI_Z!=18";
    TString selection37K = "TrackInfo.fTI_Z==19";
    TString selectronPro = "TrackInfo.fTI_Z==1";

    auto zmax = 20;
    auto hist1 = new TH2D("hist1", "top;z;x;" ,200,-200,200,200,-200,200);
    auto hist2 = new TH2D("hist2","side;z;y;" ,200,-200,200,200,-200,200);
    auto hist3 = new TH2D("hist3",";x;y;"     ,200,-200,200,200,-200,200);
    auto hist4 = new TH2D("hist4",";z;n;"     ,zmax,0,zmax,zmax,0,zmax);
    auto hist5 = new TH1D("hist5","37K;KE;"   ,200,0,100);
    auto hist6 = new TH1D("hist6","proton;KE;",200,0,100);

    tree -> Draw("TrackInfo.fTI_PositionX:TrackInfo.fTI_PositionZ>>hist1",select_not_beam,"goff");
    tree -> Draw("TrackInfo.fTI_PositionY:TrackInfo.fTI_PositionZ>>hist2",select_not_beam,"goff");
    tree -> Draw("TrackInfo.fTI_PositionY:TrackInfo.fTI_PositionX>>hist3",select_not_beam,"goff");
    tree -> Draw("TrackInfo.fTI_A-TrackInfo.fTI_Z:TrackInfo.fTI_Z>>hist4",select_not_beam,"goff");
    tree -> Draw("TrackInfo.fTI_Kinetic_Energy>>hist5",selection37K,"goff");
    tree -> Draw("TrackInfo.fTI_Kinetic_Energy>>hist6",selectronPro,"goff");

    LKDrawingGroup top("draw_simulation");
    top.AddHist(hist1);
    top.AddHist(hist2);
    top.AddHist(hist3);
    top.AddHist(hist4);
    top.AddHist(hist5);
    top.AddHist(hist6);
    top.Draw();
}
