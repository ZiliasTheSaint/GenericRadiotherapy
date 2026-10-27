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

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "UserScoreWriter.hh"

namespace
{  
  void CentralRange(G4int n, G4int& lo, G4int& hi)
  {
    if (n % 2 == 0) { lo = n/2 - 1; hi = n/2; }//e.g. n=4 (indices 0,1,2,3) => lo=1,hi=2; hi-lo+1 = 2 blocks; n=6 (0,1,2,3,4,5) => lo=2, hi=3
    else            { lo = n/2 - 1; hi = n/2 + 1; }//e.g. n=5 (indices 0,1,2,3,4) => lo=1, hi =3; hi-lo+1 = 3 blocks; n=9 (0,1,2,3,4,5,6,7,8) =>lo=3,hi=5 and 4
    if (lo < 0)     lo = 0;
    if (hi > n - 1) hi = n - 1;
  }
}

void RunAction::PrintAxisProfile(const char* axisName, G4int nBins, G4double voxelWidth,
                                       const std::vector<G4double>& sum,
                                       const std::vector<G4double>& sum2,
                                       G4int nVoxelsInBlock, G4int nEvents,
                                       G4double localConversion)
{
  for (G4int i = 0; i < nBins; i++) {
    //if (sum[i] == 0.0) continue;//contiunue loop by ignoring code below this line!
	//comment this if want all bins!!!!
    const G4double coord    = (-nBins + 1 + 2*i)*voxelWidth/2;
    const G4double meanDose = sum[i]/nVoxelsInBlock*localConversion/gray;

    const G4double err      = basicAnalyze(sum[i], sum2[i], nEvents);

	if (axisName=="x"){
		xcoordProfile.push_back(coord);
		projxProfile.push_back(meanDose);
	} else if (axisName=="y"){
		ycoordProfile.push_back(coord);
		projyProfile.push_back(meanDose);
	} else if (axisName=="z"){
		zcoordProfile.push_back(coord);
		projzProfile.push_back(meanDose);
	}

    G4cout << " " << axisName << "= " << coord/mm << " voxelDose = "
           << meanDose << " +/- " << err << " % " << G4endl;
  }
}

RunAction::RunAction()
	:totalEdepRun(0.), totalEdepRun2(0.)
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
phantom = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction(); 
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 	
	phantom->PrintDetectorParameters();

 G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;

    dose.clear();
	dose2.clear();
    EDEP.clear();
	EDEP2.clear();

	totalEdepRun=0.; totalEdepRun2=0.;
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{ 
   G4cout << "number of events = " << aRun->GetNumberOfEvent() << G4endl;
   G4int NbOfEvents = aRun->GetNumberOfEvent();
   if (NbOfEvents == 0) return;

   totalRunEnergyDeposit(NbOfEvents);
}

void RunAction::Fill(G4int index, G4double doseDeposit)
{
	//this is used for sum of doses in all voxels...has no meaning!
	dose[index] += doseDeposit;
	dose2[index] += doseDeposit*doseDeposit;
}

void RunAction::FillEnergy(G4int index, G4double energyDeposit)
{
	EDEP[index] += energyDeposit;
	EDEP2[index] += energyDeposit*energyDeposit;
}

//for better stats
void RunAction::FillEnergyTotal(G4double eventEnergy)
{
	totalEdepRun  += eventEnergy;
	totalEdepRun2 += eventEnergy*eventEnergy;
}

void RunAction::totalRunEnergyDeposit(G4int n){
	G4String resultS="";	
	
	std::map<G4int,G4double>::iterator i = dose.begin();
	std::map<G4int,G4double>::iterator end = dose.end();
	std::map<G4int,G4double>::iterator id2 = dose2.begin();

	std::map<G4int,G4double>::iterator ie = EDEP.begin();
	std::map<G4int,G4double>::iterator ie2 = EDEP2.begin();
	
	int nX = UserScoreWriter::getMeshSegment_X();
	int nY = UserScoreWriter::getMeshSegment_Y();
	int nZ = UserScoreWriter::getMeshSegment_Z();
	
	G4double voxelWidthX =2.0*UserScoreWriter::getMeshHalfSizeX()/nX;//1. *mm;
	G4double voxelWidthY =2.0*UserScoreWriter::getMeshHalfSizeY()/nY;//1. *mm;
	G4double voxelWidthZ =2.0*UserScoreWriter::getMeshHalfSizeZ()/nZ;//1. *mm;

	G4double voxelMaterialDensity = phantom->getPhantonDensity();	
	G4double voxelMass = voxelMaterialDensity*voxelWidthX*voxelWidthY*voxelWidthZ;//for renormalization
	//renormalization is needed; imagine a quanta depositing an energy E. Dose in whole phantom is E/M, but in 
	//an arbitrary voxel surrounding the deposition site is E/m which is much greater than actual phantom dose!!!
	
	G4double totalDose = 0.0;
	G4double totalDose2 = 0.0;
	G4double totalerr = 0.0;

	//just for testing
	G4double realtotalDose = 0.0;
	
	G4double realtotalerr = 0.0;
	G4double phantomMass = phantom->GetPhantomMass();
	G4double renorm = voxelMass/phantomMass;
	
	
	G4int loX, hiX, loY, hiY, loZ, hiZ;
	CentralRange(nX, loX, hiX);//n par => 2 indices; n odd =>3 indeces including central!
	CentralRange(nY, loY, hiY);
	CentralRange(nZ, loZ, hiZ);
	std::vector<G4double> profDoseX(nX, 0.0), profDoseX2(nX, 0.0);
	std::vector<G4double> profDoseY(nY, 0.0), profDoseY2(nY, 0.0);
	std::vector<G4double> profDoseZ(nZ, 0.0), profDoseZ2(nZ, 0.0);

	for(int x = 0; x < nX; x++) {
		for(int y = 0; y < nY; y++) {   
			for(int z = 0; z < nZ; z++) {				

				G4double xx = ( - nX + 1+ 2*x )* voxelWidthX/2; 
				G4double yy = ( - nY + 1+ 2*y )* voxelWidthY/2;
				G4double zz = ( - nZ + 1+ 2*z )* voxelWidthZ/2;
				
				G4int idx = x*nY*nZ +y*nZ+z;

				//find the index to be used for all maps:
				std::map<G4int, G4double>::iterator value = dose.find(idx);
				std::map<G4int, G4double>::iterator value2 = dose2.find(idx);

				std::map<G4int, G4double>::iterator value3 = EDEP.find(idx);
				std::map<G4int, G4double>::iterator value4 = EDEP2.find(idx);

				if (value != end){

					G4double ddep = (value->second);
					G4double ddep2 = (value2->second);
															
					totalDose=totalDose+ddep;//total voxel doses
					totalDose2=totalDose2+ddep2;

					G4double edep = (value3->second);
					G4double edep2 = (value4->second);

					realtotalDose = realtotalDose+edep;
					
					const G4bool cx = (x >= loX && x <= hiX);
					const G4bool cy = (y >= loY && y <= hiY);
					const G4bool cz = (z >= loZ && z <= hiZ);
					if (cx && cy) { profDoseZ[z] += ddep; profDoseZ2[z] += ddep2; }//central z profile
					if (cx && cz) { profDoseY[y] += ddep; profDoseY2[y] += ddep2; }//central y profile
					if (cy && cz) { profDoseX[x] += ddep; profDoseX2[x] += ddep2; }//central x profile
				}//if (value != end){
			}//z
		}//y
	}//x
	
	totalerr = basicAnalyze(totalDose, totalDose2, n);//not good for sum of sqares 
	//below is good we have square, not sum of squares!!!!
	realtotalerr = basicAnalyze(totalEdepRun, totalEdepRun2, n);
	
	realtotalDose=realtotalDose/phantomMass;//testing
	totalEdepRun = totalEdepRun/phantomMass;
		
	double conversion =1.0e12*renorm* 1.0/n;
	double conversion2 =1.0e12* 1.0/n;
	
	G4cout<<"Particle gun is placed at [cm] "<<
		(phantom->GetActualFocusToDetectorCenterDistance())/cm<<
		" from phantom center. It is at -z coordonate and particles are emitted toward positive z-axis."<<G4endl;
	G4cout<<"The gun generates a beam of radius [mm]: "<<
		(phantom->GetFrontalBeamRadius())/mm<<
		" at an angle [deg]: "<<(phantom->GetFrontalBeamAngle())/deg<<G4endl;	
	G4cout<<"Dose deposition (DDEP) is scored in each voxel of the mesh"<<G4endl;	
	G4cout<<"Voxel sizes [mm]: Xsize = "<<voxelWidthX/mm<<
		"; Ysize = "<<voxelWidthY/mm<<
		"; Zsize = "<<voxelWidthZ/mm<<G4endl;	

	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	
	G4cout<<"Renormalized total voxel dose [pGy/initial quantas]: "<<totalDose*conversion/gray<<" +/- "<<totalerr<<" % "<<G4endl;	
	//total dose of voxels has no meaning....the only meaning is the axis profiles!

	G4cout<<"Real phantom total dose [uGy/initial quantas]: "<<realtotalDose*conversion2/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	G4cout<<"Dose in phantom [uGy/initial quantas]: "<<totalEdepRun*conversion2/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	//confirmed from idx (realTotalDose) = from just total (totalEdepRun)
	/*
	PROOF:
	Renormalized total voxel dose [pGy/initial quantas]: 0.117142 +/- -1.23953e-021 % -> no meaning
	even with renormalization!!!!
	Real phantom total dose [uGy/initial quantas]: 0.0485184 +/- 4.24431 %//same 
	Dose in phantom [uGy/initial quantas]: 0.0485184 +/- 4.24431 %//same
	*/
	
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	
	G4cout<<"Voxel center coordinate on a particular central axis and VOXEL DOSE [pGy/initial quantas]: "<<G4endl;	
	
	G4cout<<"Local dose [mGy] on the central axes: mean over the central "
	      <<"transversal voxel block; error assumes independent voxels"<<G4endl;
	PrintAxisProfile("z", nZ, voxelWidthZ, profDoseZ, profDoseZ2, (hiX-loX+1)*(hiY-loY+1), n, conversion);
	PrintAxisProfile("y", nY, voxelWidthY, profDoseY, profDoseY2, (hiX-loX+1)*(hiZ-loZ+1), n, conversion);
	PrintAxisProfile("x", nX, voxelWidthX, profDoseX, profDoseX2, (hiY-loY+1)*(hiZ-loZ+1), n, conversion);
	
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	G4cout<<"Dose profile on central axis [pGy/initial quantas]: "<<G4endl;	
	//===========
	G4cout<<"x-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	
	for(int x = 0; x < nX; x++) {
		G4cout<<xcoordProfile[x]<<" "<<projxProfile[x]<<G4endl;		
	}
	
	G4cout<<"END_DATA"<<G4endl;

	G4cout<<"y-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;

	for(int y = 0; y < nY; y++) {
		G4cout<<ycoordProfile[y]<<" "<<projyProfile[y]<<G4endl;		
	}
	
	G4cout<<"END_DATA"<<G4endl;

	G4cout<<"z-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;

	for(int z = 0; z < nZ; z++) {
		G4cout<<zcoordProfile[z]<<" "<<projzProfile[z]<<G4endl;		
	}
	
	G4cout<<"END_DATA"<<G4endl;
}



