{
    TString lib1 = TString(gSystem->Getenv("NPLib_DIR"))+"/lib/libNPPhysics";
    TString lib2 = TString(gSystem->Getenv("NPLib_DIR"))+"/lib/libNPActar";
    TString libL = TString(gSystem->Getenv("LILAK_PATH"))+"/build/libLILAK";
    cout << lib1 << " " << (gSystem->Load(lib1)==0?"Good":"Bad") << endl;
    cout << lib2 << " " << (gSystem->Load(lib2)==0?"Good":"Bad") << endl;
    cout << libL << " " << (gSystem->Load(libL)==0?"Good":"Bad") << endl;
}
