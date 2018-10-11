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
#include "globals.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
//#include "BrachyFactory.hh"
//#include "BrachyFactoryLeipzig.hh"
//#include "BrachyFactoryIr.hh"
//#include "BrachyFactoryI.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{

 primaryMessenger = new PrimaryGeneratorMessenger(this);
 // Default source: iridium source 
 //factory = new BrachyFactoryIr();

  G4int numberParticles = 1;
 
  particleGun = new G4ParticleGun(numberParticles);
  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();  

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="gamma");//e-");//DEFAULT, change in run1.mac
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));//DEFAULT
  particleGun->SetParticleEnergy(0.662 *MeV);//DEFAULT
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
 //delete factory;
 delete particleGun;
 delete primaryMessenger;
}

void PrimaryGeneratorAction::ProcessFrontalBeam(){
  G4double FDC=Detector->GetFocusToDetectorCenterDistance();
  G4double phantomSizeZ = Detector->getPhantomSizeZ();//halfsize
  G4double height = Detector->getWorldSizeZ();//halfsize
  G4double beamRadius=Detector->GetFrontalBeamRadius();
  G4double teta =Detector->GetFrontalBeamAngle();
  if (FDC<height && FDC>phantomSizeZ)
	  height=FDC;
  Detector->SetActualFocusToDetectorCenterDistance(height);
  //overwrite:
  //height = Detector->getPhantomSizeZ();
  xin=0.0;
  yin=0.0;
  zin=-1.0 * height;//already halfsize!!!!
  uin=0.0;
  vin=0.0;
  win=cos(teta);
  
  G4double R2=0.0;
  while(true){
	  xin=G4UniformRand();
	  xin=(2.0 * xin - 1.0)*beamRadius;
	  yin=G4UniformRand();
	  yin=(2.0 * yin - 1.0)*beamRadius;
	  R2=xin*xin+yin*yin;

	  if(R2<=beamRadius*beamRadius && R2>0.0){
		  //shout at random but not good for practical purpose
		  //double cosphi=xin/sqrt(R2);
		  //double sinphi=yin/sqrt(R2);
		  //uin=sin(teta)*cosphi;
		  //vin=sin(teta)*sinphi;

		  //fixed parallel beam with a bias angle, say with respect to x-axis
		  //it is ok due to symmetry of the geometry!! 
		  vin=0.0;
		  uin=sqrt(1.0-win*win-vin*vin);

		  //now the weight
		  fweight=1.0;//here, no var reduction, source is indeed a gun!
		  break;
	  }
  }
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //factory -> CreatePrimaryGeneratorAction(anEvent);
	//if (sourceType==0)
	ProcessFrontalBeam();
	
	particleGun->SetParticlePosition(G4ThreeVector(xin,yin,zin));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(uin,vin,win));
 
    particleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::SwitchEnergy(G4String sourceChoice)
{
  /*G4int flag = 0;

  // Switch the energy spectrum of the photons delivered by the radiative source	
  if (sourceChoice == "Iodium")
    {
      flag=1;
      if (factory) delete factory;
    }
  switch(flag)
    {
    case 1:
      factory = new BrachyFactoryI;
      break;
    default:   
      factory = new BrachyFactoryIr; 
    }  */    
}
