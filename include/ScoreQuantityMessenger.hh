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

#ifndef ScoreQuantityMessenger_h
#define ScoreQuantityMessenger_h 1

#include "G4ScoreQuantityMessenger.hh"

class ScoringManager;
class G4ScoringManager;

class ScoreQuantityMessenger : public G4ScoreQuantityMessenger
{
public:
  ScoreQuantityMessenger(G4ScoringManager * SManager);
  virtual ~ScoreQuantityMessenger();
  
  void SetNewValue(G4UIcommand * command,G4String newValues);
  void MeshBinCommand(G4VScoringMesh* mesh, G4TokenVec& token);
private:
    G4ScoringManager*        fSMan;
	G4UIcommand*   qdoseDep2Cmd;
	G4UIcommand*   mBinCmd;
};
#endif