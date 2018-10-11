//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "BrachyAnalysis.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "QGSP_BIC_EMY.hh"
#include "G4PhysListFactory.hh"
#include "BrachyPhysicsList.hh"
#include "MyModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include "RunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "EventAction.hh"
#include "G4ScoringManager.hh"

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4ScoringManager.hh"
#include "UserScoreWriter.hh"
//=========
#include "ScoringManager.hh"

int main(int argc ,char ** argv)

{
  ScoringManager::GetScoringManager(); // instanciate the interactive scoring manager
	//new ScoringManager();

  G4RunManager* pRunManager = new G4RunManager;

  // Access to the Scoring Manager pointer
  G4ScoringManager* scoringManager = ScoringManager::GetScoringManager();

  // Overwrite the default output file with user-defined one 
  scoringManager->SetScoreWriter(new UserScoreWriter());

  // Initialize the physics component
  //G4VModularPhysicsList* physicsList = new QGSP_BIC_EMY;
  //physicsList->SetVerboseLevel(1);
  //pRunManager->SetUserInitialization(physicsList);
  //G4PhysListFactory factory;  //this is the new INITIALISATION!!!!!!!!!!!!!!!!!!                                   
  //pRunManager->SetUserInitialization(factory.GetReferencePhysList("QGSP_BIC_EMY");

  //pRunManager -> SetUserInitialization(new BrachyPhysicsList);
  //REMEMBER in java: WORLD SIZE>PHANTOM SIZE!!!!..external
  pRunManager->SetUserInitialization(new MyModularPhysicsList());

  // Initialize the detector component
  DetectorConstruction  *pDetectorConstruction = new  DetectorConstruction();
  pRunManager -> SetUserInitialization(pDetectorConstruction);

  // Initialize the primary particles
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction();
  pRunManager -> SetUserAction(primary);

  RunAction *pRunAction = new RunAction();
  pRunManager -> SetUserAction(pRunAction);

  pRunManager->SetUserAction(new EventAction);
  //// Initialize the Visualization component 
#ifdef G4VIS_USE
  // Visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  //Initialize G4 kernel
  pRunManager -> Initialize();

  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  if (argc == 1)   // Define UI session for interactive mode.
    { 
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      G4cout << " UI session starts ..." << G4endl;
      UImanager -> ApplyCommand("/control/execute VisualisationMacro.mac");    
      ui -> SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }  

  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete pRunManager;

  return 0;
}
