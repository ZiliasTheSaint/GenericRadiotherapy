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
  if ( edep == 0. ) return FALSE;

  G4int idx = ((G4TouchableHistory*)
	       (aStep->GetPreStepPoint()->GetTouchable()))
               ->GetReplicaNumber(indexDepth);
  G4double cubicVolume = ComputeVolume(aStep, idx);

  

  G4double density = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
    
  G4double dose    = edep / ( density * cubicVolume );
  dose *= aStep->GetPreStepPoint()->GetWeight(); 
  //dose = dose * dose;//square!!!!!!!!!!!!!!!!!!!!!!!!!!!WORKS; if commented =>we get exactly the dose!!!!!
  //Hmm..not per event but per Step!!!! Want per event!!!!
  //it is not intended for stat!!!!
  G4int  index = GetIndex(aStep);
  EvtMap->add(index,dose);  

  //G4cout<<" core voxel volume [cm3] = "<<cubicVolume/cm3
	//  <<" rho = "<<density*cm3/g
	  //<<" index = "<<index
	  //<<G4endl;

  //////
  /*const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int i = touchable->GetReplicaNumber(2);//fDepthi);
  G4int j = touchable->GetReplicaNumber(1);//fDepthj);
  G4int k = touchable->GetReplicaNumber(0);//fDepthk);
  index =  i*fNj*fNk+j*fNk+k;*/
  /////////
  //G4cout<<"RHO [g/cm3]: "<<density*cm3/g<<G4endl;//ok!
  EventAction::Fill(index,dose);//that should do the trick!
  EventAction::FillEnergy(index,edep);
  return TRUE;
}

void PSDoseDeposition2Core::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

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
  G4VSolid* solid = 0; //G4cout<<"PARAMPAMPAM!!!!!!!!!!!!!!!!!!!"<<G4endl;	
  if(physParam)
  { // for parameterized volume
	//  G4cout<<"PARAM!!!!!!!!!!!!!!!!!!!"<<G4endl;	
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
	 // G4cout<<"ORDINARY!!!!!!!!!!!!!!!!!!!"<<G4endl;	
    solid = physVol->GetLogicalVolume()->GetSolid();

	//G4cout<<"ORDINARY!!!!!!!!!!!!!!!!!!!"<<solid->GetCubicVolume()<<G4endl;	//OK!!!!!!!!!VOXEL VOLUME
  }
  
  return solid->GetCubicVolume();
}