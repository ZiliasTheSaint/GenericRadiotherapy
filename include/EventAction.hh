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
// Authors: D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
// 

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <map>

class EventActionMessenger;

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void  EndOfEventAction(const G4Event*);    
                           
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  static void Fill(G4int index, G4double energyDeposit){
	doseEventTotal[index]+= energyDeposit;//score DDEP per Event for statistics  
  };
  static void FillEnergy(G4int index, G4double energyDeposit){
	energyEventTotal[index]+= energyDeposit;//score EDEP per Event for statistics  ...used for TOTAL DOSE real (divided by TRUE total mass)
	//in dose total it is missed mass from voxels where dose dep is zero so total dose there is sum of voxel doses not TOTAL REAL DOSE in phantom
	//note that sum of voxel doses varies from mesh to mesh but REAL TOTAL DOSE per phantom is and must be a CONSTANT!!!!
	//mesh affects only errors, more mesh (dense), means more data therefore less error!!! OK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  };    
private:
   G4int     printModulo;                             
   EventActionMessenger*  eventMessenger;   
   static std::map<G4int,G4double> doseEventTotal;//Available from PSDepositionCore
   static std::map<G4int,G4double> energyEventTotal;
   void totalEventDoseDeposit();
   void totalEventEnergyDeposit();
};

#endif

    
