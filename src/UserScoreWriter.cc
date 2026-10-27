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

#include <map>
#include <fstream>

#include "UserScoreWriter.hh"
#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "G4UnitsTable.hh"

int UserScoreWriter::meshX=0;
int UserScoreWriter::meshY=0;
int UserScoreWriter::meshZ=0;

G4double UserScoreWriter::meshHalfX=0.0;
G4double UserScoreWriter::meshHalfY=0.0;
G4double UserScoreWriter::meshHalfZ=0.0;

UserScoreWriter::UserScoreWriter(): G4VScoreWriter() //call superclass constructor!
{
	// add new units for dose 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);  

  const G4double millibecqurel = 1.e-3*becquerel;
  const G4double microbecqurel = 1.e-6*becquerel;
  const G4double nanobecqurel  = 1.e-9*becquerel;  
  const G4double picobecqurel  = 1.e-12*becquerel;
  const G4double kilobecqurel = 1.e3*becquerel;
  const G4double megabecqurel = 1.e6*becquerel;
  const G4double gigabecqurel  = 1.e9*becquerel;  
  const G4double terabecqurel  = 1.e12*becquerel;

  new G4UnitDefinition("millibecqurel", "milliBq" , "Activity", millibecqurel);
  new G4UnitDefinition("microbecqurel", "microBq" , "Activity", microbecqurel);
  new G4UnitDefinition("nanobecqurel" , "nanoBq"  , "Activity", nanobecqurel);
  new G4UnitDefinition("picobecqurel" , "picoBq"  , "Activity", picobecqurel);  
  new G4UnitDefinition("kilobecqurel" , "kiloBq"  , "Activity", kilobecqurel);
  new G4UnitDefinition("megabecqurel" , "megaBq"  , "Activity", megabecqurel);
  new G4UnitDefinition("gigabecqurel" , "gigaBq"  , "Activity", gigabecqurel);
  new G4UnitDefinition("terabecqurel" , "teraBq"  , "Activity", terabecqurel);  

  iCall=0;
}

UserScoreWriter::~UserScoreWriter() 
{;}

//DEAD CODE BELLOW!!!!
void UserScoreWriter::DumpQuantityToFile(const G4String & psName, const G4String & fileName, const G4String & option) 
{
	//results are moved and displayed elsewhere!
	return;  
}

