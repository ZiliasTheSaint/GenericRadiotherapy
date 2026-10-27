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

  static void setMeshHalfSize(G4double x, G4double y, G4double z){ meshHalfX=x; meshHalfY=y; meshHalfZ=z; };
  static G4double getMeshHalfSizeX(){ return meshHalfX; };
  static G4double getMeshHalfSizeY(){ return meshHalfY; };
  static G4double getMeshHalfSizeZ(){ return meshHalfZ; };

  // store a quantity into a file
  void DumpQuantityToFile(const G4String & psName, const G4String & fileName, const G4String & option);
  
private:
	
	int iCall;
	static G4int meshX;
	static G4int meshY;
	static G4int meshZ;

	static G4double meshHalfX;
	static G4double meshHalfY;
	static G4double meshHalfZ;
};

#endif


