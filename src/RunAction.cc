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
                                       G4double localConversion
									   //added	
									   ,double& min, double& max
									   //below really not needed!
									   ,std::vector<G4double>& coordOut, std::vector<G4double>& doseOut
									   )
{
  coordOut.clear();
  doseOut.clear();
  min=100.;
  max=0.;

  for (G4int i = 0; i < nBins; i++) {
    //if (sum[i] == 0.0) continue;//contiunue loop by ignoring code below this line!
	//comment this if want all bins!!!!
    const G4double coord    = (-nBins + 1 + 2*i)*voxelWidth/2;
    const G4double meanDose = sum[i]/nVoxelsInBlock*localConversion/gray;

    const G4double err      = basicAnalyze(sum[i], sum2[i], nEvents);

	//axisName=="x" compares pointers not arrays...may work or not depending on compiler!. Better universal way:
	/*if (axisName=="x"){
		xcoordProfile.push_back(coord);
		projxProfile.push_back(meanDose);
	} else if (axisName=="y"){
		ycoordProfile.push_back(coord);
		projyProfile.push_back(meanDose);
	} else if (axisName=="z"){
		zcoordProfile.push_back(coord);
		projzProfile.push_back(meanDose);
	}*/

	coordOut.push_back(coord);
    doseOut.push_back(meanDose);

	if (err<min && err!=0.){
		min=err;
	}
	if (err>max){
		max=err;
	}

    //G4cout << " " << axisName << "= " << coord/mm << " voxelDose = "
           //<< meanDose << " +/- " << err << " % " << G4endl;
	G4cout<<coord<<" "<<meanDose<<G4endl;
  }
}

RunAction::RunAction()
	:totalEdepRun(0.), totalEdepRun2(0.)
{
  //minRelErrX = 100.; maxRellErrX = 0.; minRelErrY = 100.; maxRellErrY=0.; minRelErrZ= 100.; maxRellErrZ=0.;
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
		
	//much simpler and faster than 300 x 3 for loop:
	for (std::map<G4int,G4double>::const_iterator it = dose.begin(); it != dose.end(); ++it){
		const G4int idx = it->first;
		if (idx < 0 || idx >= nX*nY*nZ) {
			G4Exception("RunAction::totalRunEnergyDeposit()", "Etherapy002", JustWarning,
			            "Voxel index outside the mesh: check the index mapping of PSDoseDeposition2.");
			continue;
		}

		//since idx = x*nY*nZ +y*nZ+z; let  nY = 4, nZ = 5 so ny x nZ is 20 elements. Let idx = 43
		//z = 43 % 5 = 3 (idx % nZ); y = (43 / 5) % 4 = 8 % 4 = 0; x = 43 / (4 × 5) = 43 / 20 = 2
		//rebuild: idx = 2*4*5+0*5+3 = 43 corect. 43 = 2*20+3. Rest si integer division.
		//asume z=3,y=1,x=2 cu nY=4,nZ=5: idx = 2*4*5+1*5+3=48. 48/20 = 2. Restul e 1*5+3!!!!
		//restul impartirii lui idx la nZ (numarul total de elemente de pe z) ne da pozitia lui z!

		//idx e construit: Pui toti z pana la nZ. Fixeaza aceasta axa ca axa camerelor. Pozitii ocupate: nZ.
		// Axa etajelor incepe de la 0 si se termina la nY-1 avand ny etaje (fiecare cu nZ pozitii)
		//ca sa ajungi la etajul y inseamna ca ai lasat in urma y etaje complete, adica nZ*y pozitii.
		//similar axa blocurilor x. Ca sa ajungi la blocul x ai lasat in urma x*nY*nZ pozitii.
		//codificare idx = x*nY*nZ +y*nZ+z tocmai pentru a "abuza" de impartiri intregi si rest!
		//Daca vreau x din idx, observa ca idx = x*CEVA + REST, unde CEVA = nY*nZ e urias comparat cu restul de y*nZ+z
		//deci impartirea intreaga idx/(nY*nZ) ne da x. Modulul ne da y*nZ+z = mod(idx/(nY*nZ)) sau idx%(nY*nZ).
		//Similar impartirea intreaga la nZ ne da y iar modulul pe z. z = (idx%(nY*nZ)) % nZ. 
		//Matematic (idx%(nY*nZ)) % nZ e echivalent cu idx%nZ pentru ca:
		//pentru ca (y*nZ+z)%nZ = (x*nY*nZ+y*nZ+z)=z....(y*nZ)%nZ si (x*nY*nZ)%nZ sunt ambele 0!
		//similar y = (idx%(nY*nZ))/nZ ceea ce e echivalent cu y = (idx/nZ) % nY pentru ca:
		//        (y*nZ+z)/nZ = [(x*nY*nZ+y*nZ+z)/nZ] %nY=(x*nY+y)%nY = y

		const G4int x =  idx/(nY*nZ);
		const G4int y = (idx/nZ) % nY;
		const G4int z =  idx % nZ;

		const G4double ddep  = it->second;
		const G4double ddep2 = dose2.find(idx)->second;

		//edep in same for!!!!
		const std::map<G4int,G4double>::const_iterator e  = EDEP.find(idx);
		const std::map<G4int,G4double>::const_iterator e2 = EDEP2.find(idx);
		if (e != EDEP.end() && e2 != EDEP2.end()) {
			const G4double edep  = e->second;    // suma pe evenimente a energiei in voxel
			const G4double edep2 = e2->second;   // suma patratelor pe evenimente
			
			realtotalDose = realtotalDose+edep;
		}

		totalDose=totalDose+ddep;//total voxel doses
		totalDose2=totalDose2+ddep2;
		
		const G4bool cx = (x >= loX && x <= hiX);
		const G4bool cy = (y >= loY && y <= hiY);
		const G4bool cz = (z >= loZ && z <= hiZ);
		if (cx && cy) { profDoseZ[z] += ddep; profDoseZ2[z] += ddep2; }//central z profile
		if (cx && cz) { profDoseY[y] += ddep; profDoseY2[y] += ddep2; }//central y profile
		if (cy && cz) { profDoseX[x] += ddep; profDoseX2[x] += ddep2; }//central x profile
		
	}
	/*
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
	*/

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
	//still for this simple geometry renorm SHOULD WORK. This dose and 2 below should match
	//they are not...but for water they fo...seems UpdatePhantom as it is now does not work
	//it is still physically a less than mesh phantom->rest of mesh has world material, so
	//that explain why doses match only when world = phantom material!!!!!!
	G4cout<<"Phantom total dose [pGy/initial quantas] - from edep recomputed here versus filled directly from events: "<<realtotalDose*conversion2/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	G4cout<<"Dose in phantom [pGy/initial quantas]: "<<totalEdepRun*conversion2/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	//confirmed from idx (realTotalDose) = from just total (totalEdepRun)
	/*
	PROOF: and after bug fixed updating geometry in phantom construction...all 3 are the same as they should!!
	Renormalized total voxel dose [pGy/initial quantas]: 0.00470715 +/- -9.30796e-021 %
	Real phantom total dose [uGy/initial quantas]: 0.00470715 +/- 1.70689 %
	Dose in phantom [uGy/initial quantas]: 0.00470715 +/- 1.70689 %
	*/
	
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	
	G4cout<<"Voxel center coordinate on a particular central axis and VOXEL DOSE as mean over the central [pGy/initial quantas]: "<<G4endl;	
	
	//profile is a local dose -> conversion should be without renorm to phantom!!!
	//G4cout<<"Local dose [pGy/quantas] on the central axes as mean over the central: "
	      //<<"transversal voxel block; error assumes independent voxels"
	//	  <<G4endl;

	//Populate arrays directly and print results!
	G4cout<<"x-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	PrintAxisProfile("x", nX, voxelWidthX, profDoseX, profDoseX2, (hiY-loY+1)*(hiZ-loZ+1), n, conversion2, minRelErrX, maxRelErrX, xcoordProfile, projxProfile);
	G4cout<<"END_DATA"<<G4endl;
	G4cout<<"y-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	PrintAxisProfile("y", nY, voxelWidthY, profDoseY, profDoseY2, (hiX-loX+1)*(hiZ-loZ+1), n, conversion2, minRelErrY, maxRelErrY, ycoordProfile, projyProfile);
	G4cout<<"END_DATA"<<G4endl;
	G4cout<<"z-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	PrintAxisProfile("z", nZ, voxelWidthZ, profDoseZ, profDoseZ2, (hiX-loX+1)*(hiY-loY+1), n, conversion2, minRelErrZ, maxRelErrZ,  zcoordProfile, projzProfile);
	G4cout<<"END_DATA"<<G4endl;
	
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	
	G4cout<<"Voxels MIN non null and MAX relative error on X [%]: "<<minRelErrX<<", "<<maxRelErrX<<G4endl;
	G4cout<<"Voxels MIN non null and MAX relative error on Y [%]: "<<minRelErrY<<", "<<maxRelErrY<<G4endl;
	G4cout<<"Voxels MIN non null and MAX relative error on Z [%]: "<<minRelErrZ<<", "<<maxRelErrZ<<G4endl;
	
	//G4cout<<"Dose profile on central axis [pGy/initial quantas]: "<<G4endl;	
//===========THIS SECTION IS FOR SPRINGBOOT SERVER TO COLLECT DATA.....markings x-axis, BEGIN_DATA..can be moved inside print!!!!!?

	/*G4cout<<"x-axis"<<G4endl;
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
	
	G4cout<<"END_DATA"<<G4endl;*/
}



