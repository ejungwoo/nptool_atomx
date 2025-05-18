/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ATOMX simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#ifndef jw_is_debugging
#include <string>
#define jw_is_debugging std::cout<<"\033[0;32m"<<Form("+%d %s # ><><><><><><><><><>< \033[0m",__LINE__,std::string(__FILE__).c_str()) << std::endl
#define jw_cout std::cout<<"\033[0;32m"<<Form("+%d %s # \033[0m",__LINE__,std::string(__FILE__).c_str())
#endif

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ProductionCuts.hh"

// NPTool header
#include "ATOMX.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ATOMX_NS{
    // Energy and time Resolution
    const double EnergyThreshold = 0.1*MeV;
    const double ResoTime = 4.5*ns ;
    const double ResoEnergy = 1.0*MeV ;
    const double XATChamber = 1000*mm ;
    const double YATChamber = 1000*mm ;
    const double ZATChamber = 1000*mm ;

    const double Mylar_Rmax = 3.5 * cm;
    const double Mylar_Thickness = 7 * micrometer;

    const double XGasVolume = 1000. * mm;
    const double YGasVolume = 1000. * mm;
    const double ZGasVolume = 1000. * mm;

    const double XPadVolume = 800. * mm;
    const double YPadVolume = 2.   * mm;
    const double ZPadVolume = 800. * mm;

    const double XMMSVolume = 800. * mm;
    const double YMMSVolume = 220. * mm;
    const double ZMMSVolume = 800. * mm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ATOMX Specific Method
ATOMX::ATOMX(){
    m_Event = new TATOMXData();
    m_ATOMXScorer = 0;
    m_Detector = 0;

    // RGB Color + Transparency
    m_VisChamber = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.3));
    m_VisWindows = new G4VisAttributes(G4Colour(1, 0, 0, 0.25));
    m_VisGas     = new G4VisAttributes(G4Colour(0, 0.5, 0.5, 0.3));
    m_VisPads    = new G4VisAttributes(G4Colour(255, 223, 50, 0.8));
    m_VisMMS     = new G4VisAttributes(G4Colour(100, 100, 100, 0.4));
    m_VisPads -> SetForceWireframe(true);

    m_ReactionRegion = NULL;
}

ATOMX::~ATOMX(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ATOMX::AddDetector(G4ThreeVector POS, string  Shape){
    // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
    m_R.push_back(POS.mag());
    m_Theta.push_back(POS.theta());
    m_Phi.push_back(POS.phi());
    m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ATOMX::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
    m_R.push_back(R);
    m_Theta.push_back(Theta);
    m_Phi.push_back(Phi);
    m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* ATOMX::BuildDetector()
{
    G4Material* Cu = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
    G4Material* Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* Mylar = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");

    if (!m_Detector)
    {
        G4Box* sChamber    = new G4Box("ATOMX_Box", ATOMX_NS::XATChamber * 0.5, ATOMX_NS::YATChamber * 0.5, ATOMX_NS::ZATChamber * 0.5);
        G4Tubs* sWindows   = new G4Tubs("ATOMX_Windows", 0, ATOMX_NS::Mylar_Rmax, ATOMX_NS::Mylar_Thickness * 0.5, 0 * deg, 360 * deg);
        G4Box* sGas        = new G4Box("ATOMX_Gas", ATOMX_NS::XGasVolume * 0.5, ATOMX_NS::YGasVolume * 0.5, ATOMX_NS::ZGasVolume * 0.5);
        G4Box* sPad        = new G4Box("ATOMX_Pad", ATOMX_NS::XPadVolume * 0.5, ATOMX_NS::YPadVolume * 0.5, ATOMX_NS::ZPadVolume * 0.5);
        G4Box* sMMS        = new G4Box("ATOMX_MMS", ATOMX_NS::XMMSVolume * 0.5, ATOMX_NS::YMMSVolume * 0.5, ATOMX_NS::ZMMSVolume * 0.5);

        unsigned const int NumberOfGasMix = m_GasMaterial.size();

        double density = 0;
        double density_sum = 0;
        vector<G4Material*> GasComponent;
        vector<double> FractionMass;

        for (unsigned int i = 0; i < NumberOfGasMix; i++) {
            GasComponent.push_back(
                    MaterialManager::getInstance()->GetGasFromLibrary(m_GasMaterial[i], m_Pressure, m_Temperature));
        }
        for (unsigned int i = 0; i < NumberOfGasMix; i++) {
            density += ((double)m_GasFraction[i] / 100) * GasComponent[i]->GetDensity();
            density_sum += GasComponent[i]->GetDensity();
        }
        jw_cout << "density = " << density*cm3/g << endl;

        for (unsigned int i = 0; i < NumberOfGasMix; i++) {
            FractionMass.push_back(GasComponent[i]->GetDensity() / density_sum);
        }

        G4Material* GasMaterial = new G4Material("GasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);
        G4Material* DriftGasMaterial =
            new G4Material("DriftGasMix", density, NumberOfGasMix, kStateGas, m_Temperature, m_Pressure);

        for (unsigned int i = 0; i < NumberOfGasMix; i++) {
            GasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
            DriftGasMaterial->AddMaterial(GasComponent[i], FractionMass[i]);
            jw_cout << GasComponent[i] << endl;
        }

        m_Detector = new G4LogicalVolume(sChamber, GasMaterial, "logic_ATOMX_Box", 0, 0, 0);
        m_logicGas = new G4LogicalVolume(sGas, DriftGasMaterial, "logic_Gas", 0, 0, 0);
        G4LogicalVolume* logicPad = new G4LogicalVolume(sPad, Cu, "logic_Pad", 0, 0, 0);
        G4LogicalVolume* logicMMS = new G4LogicalVolume(sMMS, Al, "logic_MMS", 0, 0, 0);
        G4LogicalVolume* logicWindows = new G4LogicalVolume(sWindows, Mylar, "logic_Windows", 0, 0, 0);

        G4RotationMatrix* Rot = new G4RotationMatrix();
        new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, 0, 0)), m_logicGas, "ATOMXGas", m_Detector, false, 0);
        new G4PVPlacement(G4Transform3D(*Rot, G4ThreeVector(0, ATOMX_NS::YGasVolume * 0.5, 0)), logicPad, "ATOMXPad", m_logicGas, false, 0);

        //G4ElectricField* field = new G4UniformElectricField(G4ThreeVector(0.0, -70 * volt / cm, 0.0));
        //G4EqMagElectricField* Equation = new G4EqMagElectricField(field);
        //G4MagIntegratorStepper* Stepper = new G4ClassicalRK4(Equation, 8);

        //G4FieldManager* FieldManager = new G4FieldManager();
        //FieldManager->SetDetectorField(field);
        //m_logicGas->SetFieldManager(FieldManager, true);

        //G4MagInt_Driver* IntgrDriver = new G4MagInt_Driver(0.1 * mm, Stepper, Stepper->GetNumberOfVariables());

        //G4ChordFinder* ChordFinder = new G4ChordFinder(IntgrDriver);
        //FieldManager->SetChordFinder(ChordFinder);

        logicPad -> SetSensitiveDetector(m_ATOMXScorer);

        m_Detector  -> SetVisAttributes(m_VisChamber);
        m_logicGas  -> SetVisAttributes(m_VisGas);
        logicWindows-> SetVisAttributes(m_VisWindows);
        logicPad    -> SetVisAttributes(m_VisPads);
        logicMMS    -> SetVisAttributes(m_VisMMS);
    }

    return m_Detector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ATOMX::ReadConfiguration(NPL::InputParser parser){
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ATOMX");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl; 

    vector<string> cartSquare = {"POS", "Shape", "GasMaterial", "GasFraction", "Temperature", "Pressure"};

    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        if(blocks[i]->HasTokenList(cartSquare)){
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  ATOMX " << i+1 <<  endl;

            G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
            string Shape = blocks[i]->GetString("Shape");

            //
            vector<string> GasName = blocks[i]->GetVectorString("GasMaterial");
            vector<int> GasFraction = blocks[i]->GetVectorInt("GasFraction");
            for (unsigned int j = 0; j < GasName.size(); j++) {
                m_GasMaterial.push_back(GasName[j]);
                m_GasFraction.push_back(GasFraction[j]);
            }
            m_Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
            m_Pressure = blocks[i]->GetDouble("Pressure", "bar");
            //

            AddDetector(Pos,Shape);
        }
        else{
            cout << "ERROR: check your input file formatting " << endl;
            exit(1);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void ATOMX::ConstructDetector(G4LogicalVolume* world){
    for (unsigned short i = 0 ; i < m_R.size() ; i++) {

        G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
        G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
        G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
        G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
        // So the face of the detector is at R instead of the middle
        Det_pos+=Det_pos.unit()*ATOMX_NS::ZATChamber*0.5;
        // Building Detector reference frame
        G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
        G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
        G4double kk = -sin(m_Theta[i]);
        G4ThreeVector Y(ii,jj,kk);
        G4ThreeVector w = Det_pos.unit();
        G4ThreeVector u = w.cross(Y);
        G4ThreeVector v = w.cross(u);
        v = v.unit();
        u = u.unit();

        G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

        if(m_Shape[i] == "Square"){
            new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
                    BuildDetector(),
                    "ATOMX",world,false,i+1);
        }
    }
    if (!m_ReactionRegion) {
        G4ProductionCuts* ecut = new G4ProductionCuts();
        ecut->SetProductionCut(1000, "e-");

        m_ReactionRegion = new G4Region("NPSimulationProcess");
        m_ReactionRegion->SetProductionCuts(ecut);
        m_ReactionRegion->AddRootLogicalVolume(m_logicGas);
        m_ReactionRegion->SetUserLimits(new G4UserLimits(1.2 * mm));

        G4Region* Region_cut = new G4Region("RegionCut");
        Region_cut->SetProductionCuts(ecut);
        Region_cut->AddRootLogicalVolume(m_Detector);
    }
    G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for (unsigned int i = 0; i < size; i++) {
        mng->RemoveFastSimulationModel(m_ReactionModel[i]);
    }
    m_ReactionModel.clear();
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction", m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
    fsm = new NPS::Decay("Decay", m_ReactionRegion);
    m_ReactionModel.push_back(fsm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ATOMX::InitializeRootOutput()
{
    RootOutput *pAnalysis = RootOutput::getInstance();
    TTree *pTree = pAnalysis->GetTree();
    if(!pTree->FindBranch("ATOMX")){
        pTree->Branch("ATOMX", "TATOMXData", &m_Event) ;
    }
    pTree->SetBranchAddress("ATOMX", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ATOMX::ReadSensitive(const G4Event* ){
    m_Event->Clear();

    ///////////
    // Calorimeter scorer
    CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ATOMXScorer->GetPrimitive(0);

    unsigned int size = Scorer->GetMult(); 
    for(unsigned int i = 0 ; i < size ; i++){
        vector<unsigned int> level = Scorer->GetLevel(i); 
        double Energy = RandGauss::shoot(Scorer->GetEnergy(i),ATOMX_NS::ResoEnergy);
        if(Energy>ATOMX_NS::EnergyThreshold){
            double Time = RandGauss::shoot(Scorer->GetTime(i),ATOMX_NS::ResoTime);
            int DetectorNbr = level[0];
            m_Event->SetEnergy(DetectorNbr,Energy);
            m_Event->SetTime(DetectorNbr,Time); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void ATOMX::InitializeScorers() { 
    // This check is necessary in case the geometry is reloaded
    bool already_exist = false; 
    m_ATOMXScorer = CheckScorer("ATOMXScorer",already_exist) ;

    if(already_exist) 
        return ;

    // Otherwise the scorer is initialised
    vector<int> level; level.push_back(0);
    G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
    G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
    //and register it to the multifunctionnal detector
    m_ATOMXScorer->RegisterPrimitive(Calorimeter);
    m_ATOMXScorer->RegisterPrimitive(Interaction);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_ATOMXScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ATOMX::Construct(){
    return  (NPS::VDetector*) new ATOMX();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
    class proxy_nps_ATOMX{
        public:
            proxy_nps_ATOMX(){
                NPS::DetectorFactory::getInstance()->AddToken("ATOMX","ATOMX");
                NPS::DetectorFactory::getInstance()->AddDetector("ATOMX",ATOMX::Construct);
            }
    };

    proxy_nps_ATOMX p_nps_ATOMX;
}
