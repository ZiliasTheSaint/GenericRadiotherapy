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
#include "PSDoseDeposition2Core.hh"

#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "EventAction.hh"

PSDoseDeposition2Core::PSDoseDeposition2Core(G4String name, G4int depth)
:G4VPrimitiveScorer(name,depth),HCID(-1)
{
	SetUnit("Gy");
}

PSDoseDeposition2Core::PSDoseDeposition2Core(G4String name, const G4String& unit,
				 G4int depth)
:G4VPrimitiveScorer(name,depth),HCID(-1)
{
	SetUnit(unit);
}
PSDoseDeposition2Core::~PSDoseDeposition2Core()
{
  ;
}

G4bool PSDoseDeposition2Core::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( edep == 0. ) return false;

  const G4StepPoint* pre = aStep->GetPreStepPoint();
  const G4double weight  = pre->GetWeight();

  G4int idx = ((G4TouchableHistory*)
	       (pre->GetTouchable()))
               ->GetReplicaNumber(indexDepth);

  G4double cubicVolume = ComputeVolume(aStep, idx);
  
  G4double density = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
    
  G4double dose    = weight*edep / ( density * cubicVolume );
  
  G4int  index = GetIndex(aStep);
  EvtMap->add(index,dose);  
    
  EventAction::Fill(index,dose);
  EventAction::FillEnergy(index,weight*edep);
  return true;
}

void PSDoseDeposition2Core::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

//other methods we must override (they are virtual in superclass):
void PSDoseDeposition2Core::EndOfEvent(G4HCofThisEvent*)
{;}

void PSDoseDeposition2Core::clear()
{
  EvtMap->clear();
}

void PSDoseDeposition2Core::DrawAll()
{;}

void PSDoseDeposition2Core::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  dose deposit: " 
	   << *(itr->second)/GetUnitValue()
	   << " ["<<GetUnit() <<"]"
	   << G4endl;
  }
}

void PSDoseDeposition2Core::SetUnit(const G4String& unit)
{
	CheckAndSetUnit(unit,"Dose");
}

G4double PSDoseDeposition2Core::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0; 	
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSDoseDeposit::ComputeVolume","DetPS0004",JustWarning,ED);
    }
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume	
    solid = physVol->GetLogicalVolume()->GetSolid();
  }
  
  return solid->GetCubicVolume();
}