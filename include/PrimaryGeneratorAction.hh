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


#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
//class BrachyFactory;
class  PrimaryGeneratorMessenger;
class DetectorConstruction;
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
   PrimaryGeneratorAction();
   ~PrimaryGeneratorAction();

 public:
  void GeneratePrimaries(G4Event* anEvent);
  G4double GetParticleWeight() { return fweight;};


  void SwitchEnergy(G4String);

private:
  //BrachyFactory* factory;
  PrimaryGeneratorMessenger* primaryMessenger;
  G4ParticleGun*           particleGun;	 //pointer a to G4  class
  DetectorConstruction*    Detector;     //pointer to the geometry

  G4double fweight;//particle weight used in variance reduction
  G4double xin;G4double yin;G4double zin;//initial coordinates
  G4double uin;G4double vin;G4double win;//initial directional cosine
  void ProcessFrontalBeam();
};

#endif


