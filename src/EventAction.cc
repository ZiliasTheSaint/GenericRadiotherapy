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
// Authors: D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
// 

#include "EventAction.hh"


#include "EventActionMessenger.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include <iomanip>

std::map<G4int,G4double> EventAction::doseEventTotal;
std::map<G4int,G4double> EventAction::energyEventTotal;

EventAction::EventAction()
{  
  eventMessenger = new EventActionMessenger(this);
  printModulo = 1000;
  
}

EventAction::~EventAction()
{
  delete eventMessenger;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
  }
 
	doseEventTotal.clear();
	energyEventTotal.clear();
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  totalEventDoseDeposit();
  totalEventEnergyDeposit();
}  

void EventAction::totalEventDoseDeposit(){
	G4RunManager* runManager = G4RunManager::GetRunManager();//pointer to the static run manager
	RunAction* pointerRun = (RunAction*)(runManager->GetUserRunAction());//pointer to run action object

	std::map<G4int,G4double>::iterator i = doseEventTotal.begin();
    std::map<G4int,G4double>::iterator end = doseEventTotal.end();

	while(i!=end)
    {

       G4int index = i->first;
       G4double doseDep = i->second;
      
       //if(doseDep != 0.)
	   //{
	       pointerRun->Fill(index, doseDep);
	   //}
       
	   i++;
    }
}

void EventAction::totalEventEnergyDeposit(){
	G4RunManager* runManager = G4RunManager::GetRunManager();//pointer to the static run manager
	RunAction* pointerRun = (RunAction*)(runManager->GetUserRunAction());//pointer to run action object

	std::map<G4int,G4double>::iterator i = energyEventTotal.begin();
    std::map<G4int,G4double>::iterator end = energyEventTotal.end();

	while(i!=end)
    {

       G4int index = i->first;
       G4double edepDep = i->second;
      
      // if(edepDep != 0.)
	   //{
	       pointerRun->FillEnergy(index, edepDep);
	   //}
       
	   i++;
    }
}
