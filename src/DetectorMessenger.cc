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
//
// 

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

DetectorMessenger::DetectorMessenger( DetectorConstruction* Det): detector(Det)
{ 
  detectorDir = new G4UIdirectory("/phantom/");
  detectorDir -> SetGuidance(" phantom control.");
      
  phantomMaterialCmd = new G4UIcmdWithAString("/phantom/selectMaterial",this);
  phantomMaterialCmd -> SetGuidance("Select Material of the phantom.");
  phantomMaterialCmd -> SetParameterName("choice",false);
  phantomMaterialCmd -> AvailableForStates(G4State_Idle);
  
  mBoxSizeCmd = new G4UIcmdWith3VectorAndUnit("/other/mesh/boxSize",this);
  mBoxSizeCmd->SetGuidance("Define size of the scoring mesh.");
  mBoxSizeCmd->SetGuidance("Dx  Dy  Dz  unit");
  mBoxSizeCmd->SetParameterName("Di","Dj","Dk",false,false);
  mBoxSizeCmd->SetRange("Di>0. && Dj>0. && Dk>0.");
  mBoxSizeCmd->SetDefaultUnit("mm");

  worldSizeCmd = new G4UIcmdWithADoubleAndUnit("/world/size",this);
  worldSizeCmd->SetGuidance("Set world size");
  worldSizeCmd->SetParameterName("Size",false);
  worldSizeCmd->SetUnitCategory("Length");
  worldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fsourceTopCmd = new G4UIcmdWithADoubleAndUnit("/det/SetFocusToDetectorCenterDistance",this);
  fsourceTopCmd->SetGuidance("Define point source (or frontal beam) to detector distance.");
  fsourceTopCmd->SetParameterName("Size",false);
  fsourceTopCmd->SetRange("Size>0.");
  fsourceTopCmd->SetUnitCategory("Length");
  fsourceTopCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ffrontalBeamRadiusCmd = new G4UIcmdWithADoubleAndUnit("/det/frontalBeamRadius",this);
  ffrontalBeamRadiusCmd->SetGuidance("Define frontal beam radius.");
  ffrontalBeamRadiusCmd->SetParameterName("Size",false);
  ffrontalBeamRadiusCmd->SetRange("Size>0.");
  ffrontalBeamRadiusCmd->SetUnitCategory("Length");
  ffrontalBeamRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ffrontalBeamAngleCmd = new G4UIcmdWithADoubleAndUnit("/det/frontalBeamAngle",this);
  ffrontalBeamAngleCmd->SetGuidance("Define frontal beam angle.");
  ffrontalBeamAngleCmd->SetParameterName("Angle",false);
  ffrontalBeamAngleCmd->SetUnitCategory("Angle");
  ffrontalBeamAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 }
DetectorMessenger::~DetectorMessenger()
{
  delete phantomMaterialCmd; 
  delete mBoxSizeCmd;
  delete detectorDir;
 
  delete worldSizeCmd;//>phantom size!!!
  delete fsourceTopCmd;
	delete ffrontalBeamRadiusCmd;
	delete ffrontalBeamAngleCmd;

}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  // Change the material of the phantom
  if( command == phantomMaterialCmd )
   { detector -> SetPhantomMaterial(newValue);}

	if(command==mBoxSizeCmd) {
		//just for the command to exist for springboot server which create mac files!
		//when command will be removed from there, it will also be removed from here.
		return;
	}

	if(command==worldSizeCmd) {
		detector->setWorldSizeX(worldSizeCmd->GetNewDoubleValue(newValue));
		detector->setWorldSizeY(worldSizeCmd->GetNewDoubleValue(newValue));
		detector->setWorldSizeZ(worldSizeCmd->GetNewDoubleValue(newValue));
	}

	if( command == fsourceTopCmd ) {
    detector
      ->SetFocusToDetectorCenterDistance(fsourceTopCmd->GetNewDoubleValue(newValue));
  }

  if( command == ffrontalBeamRadiusCmd ) {
    detector
      ->SetFrontalBeamRadius(ffrontalBeamRadiusCmd->GetNewDoubleValue(newValue));
  }

  if( command == ffrontalBeamAngleCmd ) {
    detector
      ->SetFrontalBeamAngle(ffrontalBeamAngleCmd->GetNewDoubleValue(newValue));
  }
}

