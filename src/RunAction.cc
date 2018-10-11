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

struct voxelDose {
  G4double xc;//center x-coordonate
  G4double yc;//center y-coordonate
  G4double zc;//center z-coordonate
  G4double dose;//in macro filea change energyDeposit with doseDeposit
  G4double err;//percent
};

RunAction::RunAction()
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
	
	//BrachyDetectorConstruction::nbEventInRun = aRun->GetNumberOfEventToBeProcessed();
	
	//G4double rho = phantom->getPhantonDensity();
	//rho=rho*cm3/g;
	//G4cout << "Phantom density [g/cm3] = " << rho << G4endl;
 G4cout << "### Run " << aRun -> GetRunID() << " start." << G4endl;

	//doseTotal.clear();	
    dose.clear();
	dose2.clear();
    EDEP.clear();
	EDEP2.clear();
 
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
	//doseTotal[index] += doseDeposit;
	dose[index] += doseDeposit;
	dose2[index] += doseDeposit*doseDeposit;
 
	//totalPhantomEnergy += energyDeposit;
	//totalPhantomEnergy2 += energyDeposit*energyDeposit;
}

void RunAction::FillEnergy(G4int index, G4double energyDeposit)

{
	EDEP[index] += energyDeposit;
	EDEP2[index] += energyDeposit*energyDeposit;
}

void RunAction::totalRunEnergyDeposit(G4int n){
	G4String resultS="";	
	
	std::map<G4int,G4double>::iterator i = dose.begin();
	std::map<G4int,G4double>::iterator end = dose.end();
	//std::map<G4int,G4double>::iterator id = dose.begin();
	std::map<G4int,G4double>::iterator id2 = dose2.begin();

	std::map<G4int,G4double>::iterator ie = EDEP.begin();
	std::map<G4int,G4double>::iterator ie2 = EDEP2.begin();
	
	//also want initial ACTIVITY and exposureTime -1 means undefinite..until source is dead, i.e. activity<0.1 Bq!!
	/*G4double initialActivity = BrachyDetectorConstruction::initialSourceActivity;
	G4double halfLife = BrachyDetectorConstruction::sourceHalfLife;
	G4double exposureTime = BrachyDetectorConstruction::exposureTime;
	G4double zeroActivity = 0.1 * becquerel;
    if (exposureTime <0.0) {
	     //until source is exhausted 0.1Bq=A=A0 x exp(-exposureTime x ln(2)/(halfLife))
	     exposureTime = (halfLife/log(2.0))*log(initialActivity/zeroActivity);
    }
    G4double cor = (exposureTime * log(2.0) / halfLife)
						/ (1.0 - exp(-exposureTime* log(2.0) / halfLife));
    G4double meanActivityDuringExposureTime = initialActivity/cor;
    G4double NumberOfQuantaEmitted = exposureTime*meanActivityDuringExposureTime;//gamas or always present X-ray!
    //G4double NumberOfQuantaSimulated = n;//1.0 * BrachyDetectorConstruction::nbEventInRun;
	*/
	std::vector<G4double> centerZ, finalProfileZ,centerY, finalProfileY,centerX, finalProfileX;
    std::vector<voxelDose>profileZ,profileX,profileY;//, profileZmm,profileZmp,profileZpm,profileZpp;
	//std::vector<voxelDose> projX,projY,projZ;    
    //std::vector<G4double> projcenterZ, finalprojProfileZ,
		//projcenterY, finalprojProfileY,
		//projcenterX, finalprojProfileX;
    
    G4double minZ=0.0;
    G4double maxZ=0.0;
    G4double minY=0.0;
    G4double maxY=0.0;
    G4double minX=0.0;
    G4double maxX=0.0;

	int nX = UserScoreWriter::getMeshSegment_X();
	int nY = UserScoreWriter::getMeshSegment_Y();
	int nZ = UserScoreWriter::getMeshSegment_Z();
	
	G4double voxelWidthX =2.0*(phantom->getPhantomSizeX())/nX;//1. *mm;
	G4double voxelWidthY =2.0*(phantom->getPhantomSizeY())/nY;//1. *mm;
	G4double voxelWidthZ =2.0*(phantom->getPhantomSizeY())/nZ;//1. *mm;

	G4double voxelMaterialDensity = phantom->getPhantonDensity();

	//G4cout<<" voxel volume [cm3] = "<<(voxelWidthX*voxelWidthY*voxelWidthZ)/cm3
		//<<"; rho= "<<voxelMaterialDensity*cm3/g<<G4endl;

	G4double voxelMass = voxelMaterialDensity*voxelWidthX*voxelWidthY*voxelWidthZ;//for renormalization
	//renormalization is needed; imagine a quanta depositing an energy E. Dose in whole phantom is E/M, but in 
	//an arbitrary voxel surrounding the deposition site is E/m which is much greater than actual phantom dose!!!
	//renormalization, i.e. get rid of infinities, is a must, not an error but an intrinsecaly requirement.
	//it is simply how it works when dealing with small quantities such voxels. It is very similar to quantum mechanic renormalization.
	//it arise from same reasons since QM is itself a framework dealing with small quantities.
	//tested and yes, it works within resonable limits voxel total dose~real phantom dose. 
	//small differences still exists but is due to interface air(world)/tissue(phantom). !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//TESTED PUT WORLD AS BEING WATER (SAME RHO AS TISSUE) AND BOTH VOXEL TOTAL DOSE AND PHANTOM DOSE (FROM EDEP) MATCH!!! EXACTLY!!!!!!!!!!!!!!!!!!!!!

	//projection vectors
	// declare xy array
    std::vector<double> projy;std::vector<double> ycoord;
    for(int y = 0; y < nY; y++){
		G4double yy = ( - nY + 1+ 2*y )* voxelWidthY/2;
		ycoord.push_back(yy);
		projy.push_back(0.);
	}
    std::vector<double> projx;std::vector<double> xcoord;
    for(int x = 0; x < nX; x++){
		G4double xx = ( - nX + 1+ 2*x )* voxelWidthX/2;
		xcoord.push_back(xx);
		projx.push_back(0.);
	}
    std::vector<double> projz;std::vector<double> zcoord;
    for(int z = 0; z < nZ; z++){
		G4double zz = ( - nZ + 1+ 2*z )* voxelWidthZ/2;
		zcoord.push_back(zz);
		projz.push_back(0.);
	}

	std::vector<std::vector<double> > projxy;//not used
    for(int x = 0; x < nX; x++) projxy.push_back(projy);
	//
	//G4cout<<"nX = "<<nX<<" ; nY= "<<nY<<" ; nZ= "<<nZ<<G4endl;
	G4double totalDose = 0.0;
	G4double totalDose2 = 0.0;
	G4double totalerr = 0.0;//basicAnalyze(ddep, ddep, n);

	G4double realtotalDose = 0.0;
	G4double realtotalDose2 = 0.0;
	G4double realtotalerr = 0.0;//basicAnalyze(ddep, ddep, n);
	G4double phantomMass = phantom->GetPhantomMass();
	G4double renorm = voxelMass/phantomMass;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	for(int x = 0; x < nX; x++) {
		for(int y = 0; y < nY; y++) {   
			for(int z = 0; z < nZ; z++) {				

				G4double xx = ( - nX + 1+ 2*x )* voxelWidthX/2; 
				G4double yy = ( - nY + 1+ 2*y )* voxelWidthY/2;
				G4double zz = ( - nZ + 1+ 2*z )* voxelWidthZ/2;
				//if (x==0 && y==0) zcoord.push_back(zz);	
				//G4int idx = GetIndex(x, y, z);//MAPPPPPPPPPPPPPPPPPPPPPPPPPP index with value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//return x*fNMeshSegments[1]*fNMeshSegments[2] +y*fNMeshSegments[2]+z;
				G4int idx = x*nY*nZ +y*nZ+z;
				//find the index to be used for all maps:
				std::map<G4int, G4double>::iterator value = dose.find(idx);//i=doseTotal.find(idx);
				std::map<G4int, G4double>::iterator value2 = dose2.find(idx);//error

				std::map<G4int, G4double>::iterator value3 = EDEP.find(idx);//error
				std::map<G4int, G4double>::iterator value4 = EDEP2.find(idx);//error

				if (value != end){

					G4double ddep = (value->second);
					G4double ddep2 = (value2->second);
					//Note in eventAction: if(doseDep != 0.) we only have non-zero doseDep!
										
					totalDose=totalDose+ddep;//total
					totalDose2=totalDose2+ddep2;//total

					G4double edep = (value3->second);
					G4double edep2 = (value4->second);
					realtotalDose = realtotalDose+edep;///phantomMass;
					realtotalDose2 = realtotalDose2+edep2;///(phantomMass*phantomMass);
					//=======================================================================
					G4double err = basicAnalyze(ddep, ddep2, n);
					//projection profile
					projxy[x][y] += ddep;////projection on x-y plane
					projx[x]+= ddep;//projection on x
					projy[y]+= ddep;//projection on y
					projz[z]+= ddep;//projection on z
					//sum of projxyz[i] is totalDose..obviouse!

					//central axis dose
					if (
						((x==nX/2-1)&& (y==nY/2-1)) ||//((x==149)&& (y==149)) ||
						((x==nX/2-1)&& (y==nY/2)) ||//((x==149)&& (y==150)) ||
						((x==nX/2)&& (y==nY/2-1)) ||//((x==150)&& (y==149)) ||
						((x==nX/2)&& (y==nY/2))		//((x==150)&& (y==150))		
						||
						((x==nX/2-1)&& (y==nY/2+1))||
						((x==nX/2)&& (y==nY/2+1))||
						((x==nX/2+1)&& (y==nY/2+1))||
						((x==nX/2+1)&& (y==nY/2))||
						((x==nX/2+1)&& (y==nY/2-1))

						) {
				
						if (zz<minZ)minZ=zz;
						if (zz>maxZ)maxZ=zz;
						
						voxelDose vD;
						vD.xc=xx;vD.yc=yy;vD.zc=zz;
						vD.dose = ddep;
						vD.err=err;
						
						profileZ.push_back(vD);				
					}

					if (
						((x==nX/2-1)&& (z==nZ/2-1))||//((x==149)&& (z==149))||
						((x==nX/2-1)&& (z==nZ/2))||//((x==149)&& (z==150))||
						((x==nX/2)&& (z==nZ/2-1))||//((x==150)&& (z==149))||
						((x==nX/2)&& (z==nZ/2))//((x==150)&& (z==150))
						||
						((x==nX/2-1)&& (z==nZ/2+1))||
						((x==nX/2)&& (z==nZ/2+1))||
						((x==nX/2+1)&& (z==nZ/2+1))||
						((x==nX/2+1)&& (z==nZ/2))||
						((x==nX/2+1)&& (z==nZ/2-1))

						) {
			
						if (yy<minY)minY=yy;
						if (yy>maxY)maxY=yy;
						
						voxelDose vD;
						vD.xc=xx;vD.yc=yy;vD.zc=zz;
						vD.dose = ddep;
						vD.err=err;

						profileY.push_back(vD);
					//G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.ddep<<G4endl;

					}


					if (
						((y==nY/2-1)&& (z==nZ/2-1))||//((y==149)&& (z==149))||
						((y==nY/2-1)&& (z==nZ/2))||//((y==149)&& (z==150))||
						((y==nY/2)&& (z==nZ/2-1))||//((y==150)&& (z==149))||
						((y==nY/2)&& (z==nZ/2))//((y==150)&& (z==150))
						||
						((y==nY/2-1)&& (z==nZ/2+1))||
						((y==nY/2)&& (z==nZ/2+1))||
						((y==nY/2+1)&& (z==nZ/2+1))||
						((y==nY/2+1)&& (z==nZ/2))||
						((y==nY/2+1)&& (z==nZ/2-1))

						){

						if (xx<minX)minX=xx;
						if (xx>maxX)maxX=xx;
						
						voxelDose vD;
						vD.xc=xx;vD.yc=yy;vD.zc=zz;
						vD.dose = ddep;
						vD.err=err;
		    
						profileX.push_back(vD);
						//G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<
						//"; Ddep = "<<G4BestUnit(vD.dose, "Dose")<<G4endl;

					}
				}//if (value != end){
			}//z
		}//y
	}//x
	
	totalerr = basicAnalyze(totalDose, totalDose2, n);
	realtotalerr = basicAnalyze(realtotalDose, realtotalDose2, n);
	realtotalDose=realtotalDose/phantomMass;
	/*
	for(int x = 0; x < nX; x++) {
		for(int y = 0; y < fNMeshSegments[1]; y++) {

			ofile << x << "," << y << ",";
			ofile << projxy[x][y] << G4endl;

		} // y
	} // x
	*/
	

	int nmm;
	nmm= profileZ.size();
	G4double number = 1.0;//NO MEDIATION..DOSE FROM ALL 9 CELLS!!!!!!!!!!!!!!!4.0;//4 beans around center  
	G4double sum = 0.0;G4double errsum = 0.0;
	//double conversion = 1000.0*1000.0*NumberOfQuantaEmitted/n;//micro
	double conversion =1.0e12*renorm* 1.0/n;//NumberOfQuantaEmitted/n;//NumberOfQuantaSimulated;
	double conversion2 =1.0e12* 1.0/n;
	G4cout<<"Voxel dose renormalization is needed to convert voxel dose to actual phantom dose, i.e. renorm = voxelMass/phantom Mass "<<G4endl;
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
	//G4cout<<"Renormalized total voxel dose [pGy/initial quantas]: "<<totalDose*conversion/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	G4cout<<"Renormalized total voxel dose [pGy/initial quantas]: "<<totalDose*conversion/gray<<" +/- "<<totalerr<<" % "<<G4endl;	
	//G4cout<<"Real phantom total dose [uGy/initial quantas]: "<<realtotalDose*conversion2/gray<<" +/- "<<realtotalerr<<" % "<<G4endl;	
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	G4cout<<"Voxel center coordinate on a particular central axis and VOXEL DOSE [pGy/initial quantas]: "<<G4endl;	
	//G4cout<<"z-axis"<<G4endl;
	//G4cout<<"BEGIN_DATA"<<G4endl;
	for (G4double z = minZ; z<=maxZ; z=z+voxelWidthZ){//1.0){//1 mm increment
		sum = 0.0;errsum=0.0;
		for(int i = 0; i < nmm; i++) {
			voxelDose vD = profileZ[i];
			if (vD.zc==z) {
				sum = sum + vD.dose;
				errsum = errsum+vD.err;//%
			}
		}
		sum = sum /number;errsum = errsum /number;
		sum=sum*conversion/gray;//Gy->conversion
G4cout<<" z= "<<z<<" voxelDose = "<<sum<<" +/- "<<errsum<<" % "<<G4endl;
		centerZ.push_back(z);
		finalProfileZ.push_back(sum);
	}
	//G4cout<<"END_DATA"<<G4endl;

	nmm= profileY.size();
	sum = 0.0;errsum=0.0;
	//G4cout<<"y-axis"<<G4endl;
	//G4cout<<"BEGIN_DATA"<<G4endl;
	for (G4double y = minY; y<=maxY; y=y+voxelWidthY){//1.0){//1 mm increment
		sum = 0.0;errsum=0.0;
		for(int i = 0; i < nmm; i++) {
			voxelDose vD = profileY[i];
			if (vD.yc==y) {
				sum = sum + vD.dose;
				errsum = errsum+vD.err;//%
			}
		}
		sum = sum /number;errsum = errsum /number;
		sum=sum*conversion/gray;//Gy->conversion
G4cout<<" y = "<<y<<" voxelDose = "<<sum<<" +/- "<<errsum<<" % "<<G4endl;
		centerY.push_back(y);
		finalProfileY.push_back(sum);
	}
	//G4cout<<"END_DATA"<<G4endl;

	nmm= profileX.size();
	sum = 0.0;errsum=0.0;
	//G4cout<<"x-axis"<<G4endl;
	//G4cout<<"BEGIN_DATA"<<G4endl;
	for (G4double x = minX; x<=maxX; x=x+voxelWidthX){//1.0){//1 mm increment
		sum = 0.0;errsum=0.0;
		for(int i = 0; i < nmm; i++) {
			voxelDose vD = profileX[i];
			if (vD.xc==x){
				sum = sum + vD.dose;
				errsum = errsum+vD.err;//%
			}
		}
		sum = sum /number;errsum = errsum /number;
		sum=sum*conversion/gray;//Gy->conversion
G4cout<<" x= "<<x<<" voxelDose = "<<sum<<" +/- "<<errsum<<" % "<<G4endl;
		centerX.push_back(x);
		finalProfileX.push_back(sum);
	}
	//G4cout<<"END_DATA"<<G4endl;
	G4cout<<"---------------------------------------------------------------------------"<<G4endl;
	G4cout<<"Renormalized dose profile, i.e. SUM(proj_x)=SUM(proj_y)=SUM(proj_z)=total dose, [pGy/initial quantas]: "<<G4endl;	
	//===========
	G4cout<<"x-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	double dummy=0.0;double tdose = totalDose*conversion/gray;
	for(int x = 0; x < nX; x++) {
		dummy = dummy +projx[x]*conversion/gray;
		if (projx[x]!=0)
			G4cout<<xcoord[x]<<" "<<projx[x]*conversion/gray<<G4endl;
	}
	//G4cout<<" DUMMYx = "<<dummy<<"; tdose= "<<tdose<<G4endl;
	G4cout<<"END_DATA"<<G4endl;

	G4cout<<"y-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	dummy=0.0;
	for(int y = 0; y < nY; y++) {
		dummy = dummy +projy[y]*conversion/gray;
		if (projy[y]!=0)
			G4cout<<ycoord[y]<<" "<<projy[y]*conversion/gray<<G4endl;
	}
	//G4cout<<" DUMMYy = "<<dummy<<"; tdose= "<<tdose<<G4endl;
	G4cout<<"END_DATA"<<G4endl;

	G4cout<<"z-axis"<<G4endl;
	G4cout<<"BEGIN_DATA"<<G4endl;
	dummy=0.0;
	for(int z = 0; z < nZ; z++) {
		dummy = dummy +projz[z]*conversion/gray;
		if (projz[z]!=0)
			G4cout<<zcoord[z]<<" "<<projz[z]*conversion/gray<<G4endl;
	}
	//G4cout<<" DUMMYz = "<<dummy<<"; tdose= "<<tdose<<G4endl;//OK!!
	G4cout<<"END_DATA"<<G4endl;
}



