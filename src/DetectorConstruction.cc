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
#include "G4SystemOfUnits.hh"
#include "G4CSGSolid.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4TransportationManager.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "BrachyMaterial.hh"

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

G4double DetectorConstruction::phantomDensity=1.0;//solved

DetectorConstruction::DetectorConstruction()
: detectorChoice(0),
  World(0), WorldLog(0), WorldPhys(0),
  Phantom(0), PhantomLog(0), PhantomPhys(0),
  phantomAbsorberMaterial(0)
{
  // Define half size of the phantom along the x, y, z axis
  phantomSizeX = 5.*cm ;
  phantomSizeY = 5.*cm;
  phantomSizeZ = 5.*cm;
  
  // Define the sizes of the World volume containing the phantom
  worldSizeX = 4.0*m;
  worldSizeY = 4.0*m;
  worldSizeZ = 4.0*m;

  fFrontalBeam_radius=0.1*cm;//source is parallel frontal beam
  fFrontalBeam_angle=0.0;//frontal beam angle to world z-axis
  fFocusToDetectorCenterDistance = 0.5*worldSizeZ;

  phantomMass = 1.0;//init

  // Define the messenger of the Detector component
  // It is possible to modify geometrical parameters through UI
  detectorMessenger = new DetectorMessenger(this);
  
  pMaterial = new BrachyMaterial();
}

DetectorConstruction::~DetectorConstruction()
{ 
  delete pMaterial;
  delete detectorMessenger;
}

void DetectorConstruction::UpdatePhantomSize(G4double hx, G4double hy, G4double hz)
{
	//fix===============================
	//see hh G4Box*              Phantom;  //pointer to solid phantom
	if (!Phantom) return;                       // geometria nu e construita inca
	if (hx > worldSizeX || hy > worldSizeY || hz > worldSizeZ)
    G4Exception("DetectorConstruction::UpdatePhantomSize()", "Etherapy001", FatalException,
                "The phantom (mesh) does not fit in the world volume.");
	//=======================================

    // update to mesh size
    phantomSizeX = hx;
    phantomSizeY = hy;
    phantomSizeZ = hz;

	//Fix=========================
	Phantom->SetXHalfLength(hx);                // solidul real
	Phantom->SetYHalfLength(hy);
	Phantom->SetZHalfLength(hz);
	phantomMass = phantomDensity * 8.0*hx*hy*hz;//new mass. 2^3=8 because half size!!!!
	//=============================

    // Tell Geant4 geometry was changed
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4cout << "--> Phantom geometry updated!!!" << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  pMaterial -> DefineMaterials();
 
  // Model the phantom (water box)  
  ConstructPhantom();

  return WorldPhys;
}

void DetectorConstruction::ConstructPhantom()
{
  // Model the water phantom 
  
  // Define the light blue color
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* air = pMaterial -> GetMat("Air") ;
  G4Material* water = pMaterial -> GetMat("Water");
  phantomDensity = water->GetDensity();

  // World volume
  World = new G4Box("World",worldSizeX,worldSizeY,worldSizeZ);
  WorldLog = new G4LogicalVolume(World,air,"WorldLog",0,0,0);      //total dose != dose='real'????
  //WorldLog = new G4LogicalVolume(World,water,"WorldLog",0,0,0);  //total dose=dose='real'??????
  /*
  
  UpdatePhantomSize() nu redimensioneaza fantomul.
  G4Box-ul Phantom construit in ConstructPhantom() nu e atins nicaieri, iar Construct() nu mai ruleaza 
  dupa Initialize(). GeometryHasBeenModified() doar re-inchide si re-optimizeaza geometria existenta.
  So: Geometria reala ramane cubul de 10 cm, iar PrintDetectorParameters() afiseazã 30 cm, pentru ca
  tipareste variabilele, nu solidul. De aceea pare cã „mesh = fantom“.
	phantomMass ramane 1 kg (cubul de 10 cm), desi mesh-ul acoperã 27 kg.
	Voxelii din afara fantomului sunt aer, cu doza ×~830 fata de apa.

	Asta raspunde la intrebarea din ultimul commit (de ce dozele sunt la fel doar cand lumea are acelasi
	material cu fantomu“). Totalul voxel e SUM E_v x (RHO_fantom/RHO_v) / M, deci coincide cu 
	SUM E / M doar daca RHO_v = RHO_fantom peste tot. Comentariul din RunAction.cc arata exact acest
	efect: 0,117 fata de 0,0485.!!!!!!!!!!
  
  */
  WorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                "WorldPhys",WorldLog,0,false,0);

  // Water Box
  Phantom = new G4Box("Phantom",phantomSizeX,phantomSizeY,
                       phantomSizeZ);
  // Logical volume
  PhantomLog = new G4LogicalVolume(Phantom,water,"PhantomLog",0,0,0);

  // Physical volume
  PhantomPhys = new G4PVPlacement(0,G4ThreeVector(), // Position: rotation and translation
                  "PhantomPhys", // Name
				  PhantomLog, // Associated logical volume
                  WorldPhys, // Mother volume
				  false,0); 
  
  G4double phantomVolume = PhantomLog->GetSolid()->GetCubicVolume();
  phantomMass = phantomDensity * phantomVolume;
  //WorldLog -> SetVisAttributes (G4VisAttributes::Invisible);

  // Visualization attributes of the phantom
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(lblue);
  simpleBoxVisAtt -> SetVisibility(true);
  simpleBoxVisAtt -> SetForceWireframe(true);
  PhantomLog -> SetVisAttributes(simpleBoxVisAtt);
}


void DetectorConstruction::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the phantom is a box whose size is (X,Y,Z): " << G4endl
         << 2.0*Phantom->GetXHalfLength()/cm//real dimension to not lie if not catch bug when not updated geometry!!phantomSizeX *2./cm
         << " cm X "
         << 2.0*Phantom->GetYHalfLength()/cm//phantomSizeY *2./cm
         << " cm X "
         << 2.0*Phantom->GetYHalfLength()/cm//phantomSizeZ *2./cm
         << " cm" << G4endl
         << "The phantom is made of "
         << phantomAbsorberMaterial -> GetName() <<G4endl
		 << "Phantom density [g/cm3] is: "<<phantomDensity*cm3/g<<G4endl         
         << "-------------------------------------------------------------------------"
         << G4endl;
}

void DetectorConstruction::SetPhantomMaterial(G4String materialChoice)
{
  // Search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
  {
    phantomAbsorberMaterial = pttoMaterial;
	phantomDensity = pttoMaterial->GetDensity();
    PhantomLog -> SetMaterial(pttoMaterial); 

	G4double phantomVolume = PhantomLog->GetSolid()->GetCubicVolume();
    phantomMass = phantomDensity * phantomVolume;
  } 
  else
    G4cout << "WARNING: material '" << materialChoice
           << "' not available!" << G4endl;            
}
