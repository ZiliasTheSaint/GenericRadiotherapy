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

#ifndef UserScoreWriter_h
#define UserScoreWriter_h 1

#include "globals.hh"
#include "G4VScoreWriter.hh"

class UserScoreWriter : public G4VScoreWriter {

public:
  UserScoreWriter();
  virtual ~UserScoreWriter();

public:
  static G4int getMeshSegment_X(){
	  return meshX;
  };
  static G4int getMeshSegment_Y(){
	  return meshY;
  };
  static G4int getMeshSegment_Z(){
	  return meshZ;
  };
  static void setMeshSegment_X(G4int val){
	  meshX = val;
  };
  static void setMeshSegment_Y(G4int val){
	  meshY = val;
  };
  static void setMeshSegment_Z(G4int val){
	  meshZ = val;
  };

  // store a quantity into a file
  void DumpQuantityToFile(const G4String & psName, const G4String & fileName, const G4String & option);
  
  G4double basicAnalyze(G4double x, G4double x2, G4int n){
	  G4double temp=x;G4double temp2=x2;G4int nevents=n;
	  temp=temp/nevents;temp2=temp2/nevents;
	  temp2 = temp2 - temp*temp;
	  if (nevents>1) temp2=temp2/(nevents-1.0);
	  if (temp2 >0.)
	   temp2 = sqrt(temp2); 
	  //else temp2 = 99.99;//never!
	  //=========percent
	  if (temp!=0.0){
	   temp2 = std::min(100.0*temp2/temp,99.9);
      } //else temp2 = 99.9;//no score means no score not necessarly an error!
	  return temp2;
  }

private:
	//BrachyDetectorConstruction*    phantom;     //pointer to the geometry
	int iCall;
	static G4int meshX;
	static G4int meshY;
	static G4int meshZ;
};

#endif


