void draw_summary()
{
    auto graph = new TGraph();
    TSystemFile *sysFile = nullptr;
    TIter nextFile(TSystemDirectory("search","data").GetListOfFiles());
    while ((sysFile=(TSystemFile*)nextFile()))
    {
        TString fileName = sysFile -> GetName();
        if (fileName.Index("efficiency")==0)
        {
            auto i1 = fileName.Index("_z");
            auto i2 = fileName.Index(".ana");
            TString z_value = fileName(i1+2,i2-i1-4);
            cout << fileName << endl;
            ifstream file(Form("data/%s",fileName.Data()));
            double eff_value;
            file >> eff_value;
            graph -> SetPoint(graph->GetN(),z_value.Atof(),eff_value);
        }
    }
    graph -> Sort();
    graph -> Print();

    auto hist = new TH2D("hist",";z (mm);geometrical efficiency",100,-160,160,100,0.18,0.55);
    hist -> SetStats(0);
    hist -> Draw();
    graph -> Draw("same*l");

    gPad -> SetGridx();
    gPad -> SetGridy();
    gPad -> SaveAs("figures/efficiency.png");
    gPad -> SaveAs("figures/efficiency.pdf");
}
