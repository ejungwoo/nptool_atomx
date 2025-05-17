{
    for (TString path : {
            TString(gSystem->Getenv("NPLib_DIR"))+"/lib/libNPPhysics",
            TString(gSystem->Getenv("NPLib_DIR"))+"/lib/libNPActar"
            }) {
        cout << path << " " << gSystem -> Load(path) << endl;
    }

    TString libString = TString(gSystem->Getenv("LILAK_PATH"))+"/build/libLILAK";
    int loadv = gSystem -> Load(libString);
    if      (loadv== 0) cout << "LILAK" << endl;
    else if (loadv== 1) cout << "LILAK (already loaded)" << endl;
    else if (loadv==-1) cout << "ERROR: Cannot Load LILAK from " << libString << endl;
    else if (loadv==-2) cout << "ERROR: LILAK version miss match!" << libString << endl;
}
