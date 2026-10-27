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
// Authors: D. Fulea, National Institute of Public Health Bucharest, Cluj-Napoca Regional Center, Romania
//
#include "ScoreQuantityMessenger.hh"
#include "ScoringManager.hh"
#include "G4ScoringManager.hh"
#include "PSDoseDeposition2.hh"
#include "UserScoreWriter.hh"

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

ScoreQuantityMessenger::ScoreQuantityMessenger(G4ScoringManager* SManager)
:G4ScoreQuantityMessenger(SManager),fSMan(SManager)//call superclass constructor and set fsMan = SManager
{
   G4UIparameter* param;

  qdoseDep2Cmd = new G4UIcommand("/score/quantity/doseDeposit2",this);
  qdoseDep2Cmd->SetGuidance("Custom dose deposit scorer.");
  qdoseDep2Cmd->
  SetGuidance("[usage] /score/quantiy/doseDeposit2 qname unit");
  qdoseDep2Cmd->SetGuidance("  qname  :(String) scorer name");
  qdoseDep2Cmd->SetGuidance("  unit   :(String) unit");
  param = new G4UIparameter("qname",'s',false);
  qdoseDep2Cmd->SetParameter(param);
  param = new G4UIparameter("unit",'s',true);
  param->SetDefaultValue("Gy");
  qdoseDep2Cmd->SetParameter(param);

  mBinCmd = new G4UIcommand("/other/mesh/nBin",this);
  mBinCmd->SetGuidance("Define segments of the scoring mesh.");
  param = new G4UIparameter("Ni",'i',false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Ni>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nj",'i',false);
  param->SetDefaultValue("1");
  param->SetParameterRange("Nj>0");
  mBinCmd->SetParameter(param);
  param = new G4UIparameter("Nk",'i',false);
  param->SetDefaultValue("1");
  mBinCmd->SetParameter(param);
  param->SetParameterRange("Nk>0");
 
}

ScoreQuantityMessenger::~ScoreQuantityMessenger()
{
  delete    qdoseDep2Cmd;
  delete    mBinCmd;
}

void ScoreQuantityMessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{	
	G4VScoringMesh* mesh = fSMan->GetCurrentMesh();
      if(!mesh)
      {
        G4cerr << "ERROR : No mesh is currently open. Open/create a mesh first. Command ignored." << G4endl;
        return;
      }
      // Tokens
      G4TokenVec token;
      FillTokenVec(newVal,token);
    
      // Commands for Current Mesh
	  if(command== qdoseDep2Cmd) {
		if ( CheckMeshPS(mesh,token[0]) ){
           if( mesh->GetShape()==boxMesh ) {
			  
			    PSDoseDeposition2* ps = new PSDoseDeposition2(token[0]);
				ps->SetUnit(token[1]);
				mesh->SetPrimitiveScorer(ps);

				G4ThreeVector msize = mesh->GetSize();
				UserScoreWriter::setMeshHalfSize(msize.x(), msize.y(), msize.z());

				G4int nSeg[3];
				mesh->GetNumberOfSegments(nSeg);
								
				UserScoreWriter::setMeshSegment_X(nSeg[0]);
				UserScoreWriter::setMeshSegment_Y(nSeg[1]);
				UserScoreWriter::setMeshSegment_Z(nSeg[2]);

				G4RunManager* runManager = G4RunManager::GetRunManager();
				DetectorConstruction* detector = (DetectorConstruction*)(runManager->GetUserDetectorConstruction());

				if(detector) {
					// phantom must match the mesh
					detector->UpdatePhantomSize(msize.x(), msize.y(), msize.z());
				 }
            } 
			else if( mesh->GetShape()==cylinderMesh ) {
			  //dead code left for reference
			  //OUR MESH IS FIXED SIZE -BOX- example of cylinder for future:

              //G4PSDoseDepositForCylinder3D* ps = 
			  //new G4PSDoseDepositForCylinder3D(token[0]);
              //ps->SetUnit(token[1]);
	          //G4ThreeVector msize = mesh->GetSize(); // gevin in R Z N/A
              //ps->SetCylinderSize(msize[0],msize[1]); // given in dr dz
              //G4int nSeg[3];
              //mesh->GetNumberOfSegments(nSeg);
              //ps->SetNumberOfSegments(nSeg);
	          //mesh->SetPrimitiveScorer(ps);
            }
		}
      } else if(command==mBinCmd) {
	      MeshBinCommand(mesh,token);
	  } 
	  
}

void ScoreQuantityMessenger::MeshBinCommand(G4VScoringMesh* mesh, G4TokenVec& token){
    //just for the command to exist. Macro files are auto generated by springboot server and
	//the /other/mesh/nBin command is written.
}
