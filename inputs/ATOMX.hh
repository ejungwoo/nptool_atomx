#ifndef ATOMX_h
#define ATOMX_h 1
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

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4UserLimits.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4VFastSimulationModel.hh"
#include "G4FastSimulationManager.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "Decay.hh"
#include "TATOMXData.h"
#include "NPSVDetector.hh"
#include "NPInputParser.h"
#include "BeamReaction.hh"

class ATOMX : public NPS::VDetector{
    ////////////////////////////////////////////////////
    /////// Default Constructor and Destructor /////////
    ////////////////////////////////////////////////////
    public:
        ATOMX() ;
        virtual ~ATOMX() ;

        ////////////////////////////////////////////////////
        /////// Specific Function of this Class ///////////
        ////////////////////////////////////////////////////
    public:
        // Cartesian
        void AddDetector(G4ThreeVector POS, string Shape);
        // Spherical
        void AddDetector(double R,double Theta,double Phi,string Shape);  


        G4LogicalVolume* BuildDetector();

    private:
        G4LogicalVolume* m_Detector;
        G4LogicalVolume* m_logicGas;

        ////////////////////////////////////////////////////
        //////  Inherite from NPS::VDetector class /////////
        ////////////////////////////////////////////////////
    public:
        // Read stream at Configfile to pick-up parameters of detector (Position,...)
        // Called in DetecorConstruction::ReadDetextorConfiguration Method
        void ReadConfiguration(NPL::InputParser) ;

        // Construct detector and inialise sensitive part.
        // Called After DetecorConstruction::AddDetector Method
        void ConstructDetector(G4LogicalVolume* world) ;

        // Add Detector branch to the EventTree.
        // Called After DetecorConstruction::AddDetector Method
        void InitializeRootOutput() ;

        // Read sensitive part and fill the Root tree.
        // Called at in the EventAction::EndOfEventAvtion
        void ReadSensitive(const G4Event* event) ;

    public:   // Scorer
        //   Initialize all Scorer used by the MUST2Array
        void InitializeScorers() ;

        //   Associated Scorer
        G4MultiFunctionalDetector* m_ATOMXScorer ;
        ////////////////////////////////////////////////////
        ///////////Event class to store Data////////////////
        ////////////////////////////////////////////////////
    private:
        TATOMXData* m_Event ;

        ////////////////////////////////////////////////////
        ///////////////Private intern Data//////////////////
        ////////////////////////////////////////////////////
    private: // Geometry
        // Detector Coordinate 
        vector<double>  m_R; 
        vector<double>  m_Theta;
        vector<double>  m_Phi; 

        //   Shape type
        vector<string> m_Shape ;

        // token
        vector<string> m_GasMaterial;
        vector<int> m_GasFraction;
        double m_Pressure; // bar
        double m_Temperature; // kelvin

        // Visualisation Attribute
        G4VisAttributes* m_VisChamber;
        G4VisAttributes* m_VisWindows;
        G4VisAttributes* m_VisMMS;
        G4VisAttributes* m_VisGas;
        G4VisAttributes* m_VisPads;

    private:
        // Region were reaction can occure:
        G4Region* m_ReactionRegion;
        vector<G4VFastSimulationModel*> m_ReactionModel;

    public:
        static NPS::VDetector* Construct();
};
#endif
