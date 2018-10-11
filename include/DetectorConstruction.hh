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

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4Material;
class G4Box;
class G4Colour;
class G4VPhysicalVolume;
class G4VPhysicalVolume;
class BrachyMaterial;
//class BrachyFactory;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume*   Construct();  
  void SwitchBrachytherapicSeed(); //Change radiactive source through GUI
  void SelectBrachytherapicSeed(G4String val);

  void ConstructPhantom(); 
  void PrintDetectorParameters(); 
  void SetPhantomMaterial(G4String); 
  G4double getPhantonDensity(){return phantomDensity;};
  static G4double phantomDensity;
  //static G4double initialSourceActivity;
  //static G4double sourceHalfLife;
  //static G4double exposureTime;
  //static G4int nbEventInRun;
  
  void setPhantomSizeX(G4double val){
	  phantomSizeX = val;//in fact halfSize
  };
  void setPhantomSizeY(G4double val){
	  phantomSizeY = val;//in fact halfSize
  };
  void setPhantomSizeZ(G4double val){
	  phantomSizeZ = val;//in fact halfSize
  };
  G4double getPhantomSizeX(){
	  return phantomSizeX;//in fact halfSize
  };
  G4double getPhantomSizeY(){
	  return phantomSizeY;//in fact halfSize
  };
  G4double getPhantomSizeZ(){
	  return phantomSizeZ;//in fact halfSize
  };
  void setWorldSizeX(G4double val){
	  worldSizeX = val;//in fact halfSize
  };
  void setWorldSizeY(G4double val){
	  worldSizeY = val;//in fact halfSize
  };
  void setWorldSizeZ(G4double val){
	  worldSizeZ = val;//in fact halfSize
  };
  G4double getWorldSizeX(){
	  return worldSizeX;//in fact halfSize
  };
  G4double getWorldSizeY(){
	  return worldSizeY;//in fact halfSize
  };
  G4double getWorldSizeZ(){
	  return worldSizeZ;//in fact halfSize
  };
  G4double GetActualFocusToDetectorCenterDistance()           {return fFocusToDetectorCenterRealDistance;};
  G4double GetFocusToDetectorCenterDistance()           {return fFocusToDetectorCenterDistance;};
  G4double GetFrontalBeamRadius()           {return fFrontalBeam_radius;};
  G4double GetFrontalBeamAngle()           {return fFrontalBeam_angle;};
  void SetFrontalBeamRadius(G4double val){fFrontalBeam_radius=val;};
  void SetFrontalBeamAngle(G4double val){fFrontalBeam_angle=val;};
  void SetFocusToDetectorCenterDistance(G4double val){fFocusToDetectorCenterDistance=val;};
  void SetActualFocusToDetectorCenterDistance(G4double val){fFocusToDetectorCenterRealDistance=val;};

  G4double GetPhantomMass()           {return phantomMass;};
private:
  
  G4int detectorChoice; //Select brachytherapic seed
  //BrachyFactory* factory;

  // World ...
  G4Box*             World;        //pointer to the solid World 
  G4LogicalVolume*   WorldLog;     //pointer to the logical World
  G4VPhysicalVolume* WorldPhys;    //pointer to the physical World

  // Phantom ... 
  G4Box*              Phantom;  //pointer to solid phantom
  G4LogicalVolume*    PhantomLog; //pointer to logic phantom
  G4VPhysicalVolume*  PhantomPhys; //pointer to physical phantom
  G4Material*         phantomAbsorberMaterial;
 
  G4double phantomSizeX; //Phantom XSize
  G4double phantomSizeY; //Phantom YSize
  G4double phantomSizeZ; //Phantom ZSize  
  G4double worldSizeX ; //World XSize
  G4double worldSizeY ; //World YSize
  G4double worldSizeZ ; //World XSize
  DetectorMessenger* detectorMessenger; 
  BrachyMaterial* pMaterial;   

  G4double fFrontalBeam_radius;
  G4double fFrontalBeam_angle;
  G4double fFocusToDetectorCenterDistance;
  G4double fFocusToDetectorCenterRealDistance;
  G4double phantomMass;
  //G4int     printModulo;  
};

#endif
