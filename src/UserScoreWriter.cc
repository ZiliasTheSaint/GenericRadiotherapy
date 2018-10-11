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
#include "BrachyAnalysis.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

#include "G4UnitsTable.hh"
//#include "G4ios.hh"

int UserScoreWriter::meshX=0;
int UserScoreWriter::meshY=0;
int UserScoreWriter::meshZ=0;

struct voxelDose {
  G4double xc;//center x-coordonate
  G4double yc;//center y-coordonate
  G4double zc;//center z-coordonate
  //G4double edep;
  G4double dose;//in macro filea change energyDeposit with doseDeposit
};

// The default output is
// voxelX, voxelY, voxelZ, edep
// The BrachyUserScoreWriter allows to change the format of the output file.
// in the specific case:
// xx (mm)  yy(mm) zz(mm) edep(keV)

// The same information is stored in a ntuple, in the 
// brachytherapy.root file

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

  //meshX=fNMeshSegments[0];
  //meshY=fNMeshSegments[1];
  //meshZ=fNMeshSegments[2];
//phantom = (BrachyDetectorConstruction*)
  //           G4RunManager::GetRunManager()->GetUserDetectorConstruction(); 
}

UserScoreWriter::~UserScoreWriter() 
{;}

void UserScoreWriter::DumpQuantityToFile(const G4String & psName, const G4String & fileName, const G4String & option) 
{
	return;
  if(verboseLevel > 0) {
    G4cout << "BrachyUserScorer-defined DumpQuantityToFile() method is invoked."
           << G4endl; }

  // change the option string into lowercase to the case-insensitive.
  G4String opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

  // confirm the option
  if(opt.size() == 0) opt = "csv";

  // open the file
  std::ofstream ofile(fileName);//std::ofstream ofile2(fileName+"_Z");
  if(!ofile) {
    G4cerr << "ERROR : DumpToFile : File open error -> "
	   << fileName << G4endl;
    return;
  }
  ofile << "# mesh name: " << fScoringMesh->GetWorldName() << G4endl;//ofile2 << "# mesh name: " << fScoringMesh->GetWorldName() << G4endl;
  
  //call
  iCall=iCall+1;//ifiCall=1 dose, if iCall=2 we have dose2
  // retrieve the map
  MeshScoreMap fSMap = fScoringMesh -> GetScoreMap();
  
  MeshScoreMap::const_iterator msMapItr = fSMap.find(psName);
  if(msMapItr == fSMap.end()) {
    G4cerr << "ERROR : DumpToFile : Unknown quantity, \""
	   << psName << "\"." << G4endl;
    return;
  }
  std::map<G4int, G4double*> * score = msMapItr -> second-> GetMap();
 
  ofile << "# primitive scorer name: " << msMapItr -> first << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //
  // Open a ROOT output file
  //

  analysisManager -> OpenFile("brachytherapy");
   
  //
  // Creating ntuple
  //
  analysisManager -> CreateNtuple("EnergyDeposition", "Edep(keV) in the phantom");
  analysisManager -> CreateNtupleDColumn("xx");
  analysisManager -> CreateNtupleDColumn("yy");
  analysisManager -> CreateNtupleDColumn("zz");
  analysisManager -> CreateNtupleDColumn("edep");
  analysisManager -> FinishNtuple();
  
 //
 // Write quantity in the ASCII output file and in brachytherapy.root
 //

  // declare xy array
  //std::vector<double> projy;
  //for(int y = 0; y < fNMeshSegments[1]; y++) projy.push_back(0.);
  //std::vector<std::vector<double> > projxy;
  //for(int x = 0; x < fNMeshSegments[0]; x++) projxy.push_back(projy);

  //G4double rho = phantom->getPhantonDensity();
/*  G4double rho = BrachyDetectorConstruction::phantomDensity;  
  G4double voxelSize = 1. *mm;
  G4double voxelVolume = voxelSize*voxelSize*voxelSize;
  G4double voxelMass =rho*voxelVolume;
  //G4cout<<"VoxelMass [mg]: "<<voxelMass*1000.0/g<<G4endl;

  //also want initial ACTIVITY and exposureTime -1 means undefinite..until source is dead, i.e. activity<0.1 Bq!!
  G4double initialActivity = BrachyDetectorConstruction::initialSourceActivity;
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
  G4double NumberOfQuantaSimulated = 1.0 * BrachyDetectorConstruction::nbEventInRun;
  //G4cout<<"MeanAct: "<<G4BestUnit(meanActivityDuringExposureTime, "Activity")<<G4endl;
  //G4cout<<"Nquantas: "<<NumberOfQuantaEmmitted<<G4endl;
  

  std::vector<G4double> centerZ, finalProfileZ,centerY, finalProfileY,centerX, finalProfileX;
  std::vector<voxelDose>profileZ,profileX,profileY;//, profileZmm,profileZmp,profileZpm,profileZpp;
  G4double minZ=0.0;
  G4double maxZ=0.0;
  G4double minY=0.0;
  G4double maxY=0.0;
  G4double minX=0.0;
  G4double maxX=0.0;

  */
 ofile << std::setprecision(16); // for double value with 8 bytes
							//ofile2 << std::setprecision(16);
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {   
       for(int z = 0; z < fNMeshSegments[2]; z++) {

       G4int numberOfVoxel = fNMeshSegments[0];

       // If the voxel width is changed in the macro file, 
       // the voxel width variable must be updated
       G4double voxelWidth = 1. *mm;
       //

       G4double xx = ( - numberOfVoxel + 1+ 2*x )* voxelWidth/2; 
       G4double yy = ( - numberOfVoxel + 1+ 2*y )* voxelWidth/2;
       G4double zz = ( - numberOfVoxel + 1+ 2*z )* voxelWidth/2;

        G4int idx = GetIndex(x, y, z);//MAPPPPPPPPPPPPPPPPPPPPPPPPPP index with value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	std::map<G4int, G4double*>::iterator value = score -> find(idx);
//=======================================================================
/*	if ((x==149)&& (y==149)) {
		if (value != score -> end()) {
			//ofile2 << zz <<"  " <<*(value->second)/keV << G4endl;
			//ofile2 << xx << "  " << yy << "  " << zz <<"  " <<*(value->second)/keV << G4endl;
			if (zz<minZ)minZ=zz;
			if (zz>maxZ)maxZ=zz;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
			vD.xc=xx;vD.yc=yy;vD.zc=zz;//vD.edep=edep;
			vD.dose = edep;//edep/voxelMass;

		    //profileZmm.push_back(vD);
			profileZ.push_back(vD);
		}
	}
if ((x==149)&& (y==150)) {
		if (value != score -> end()) {

			if (zz<minZ)minZ=zz;
			if (zz>maxZ)maxZ=zz;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
		    vD.xc=xx;vD.yc=yy;vD.zc=zz;//vD.edep=edep;
			vD.dose = edep;//edep/voxelMass;
		    //profileZmp.push_back(vD);
			profileZ.push_back(vD);
		}
	}
if ((x==150)&& (y==149)) {
		if (value != score -> end()) {
			
			if (zz<minZ)minZ=zz;
			if (zz>maxZ)maxZ=zz;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
		    vD.xc=xx;vD.yc=yy;vD.zc=zz;//vD.edep=edep;
			vD.dose = edep;//edep/voxelMass;
		    //profileZpm.push_back(vD);
			profileZ.push_back(vD);
		}
	}
if ((x==150)&& (y==150)) {
		if (value != score -> end()) {
			
			if (zz<minZ)minZ=zz;
			if (zz>maxZ)maxZ=zz;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
		    vD.xc=xx;vD.yc=yy;vD.zc=zz;
			//vD.edep=edep;
			vD.dose = edep;///voxelMass;
		    //profileZpp.push_back(vD);
			profileZ.push_back(vD);
		}
	}

if (
	((x==149)&& (z==149))||
	((x==149)&& (z==150))||
	((x==150)&& (z==149))||
	((x==150)&& (z==150))
	)
   {
		if (value != score -> end()) {
			//ofile2 << zz <<"  " <<*(value->second)/keV << G4endl;
			//ofile2 << xx << "  " << yy << "  " << zz <<"  " <<*(value->second)/keV << G4endl;
			if (yy<minY)minY=yy;
			if (yy>maxY)maxY=yy;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
		    vD.xc=xx;vD.yc=yy;vD.zc=zz;//vD.edep=edep;
			vD.dose = edep;//edep/voxelMass;
		    //profileZmm.push_back(vD);
			profileY.push_back(vD);
//G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
		}
	}

if (
	((y==149)&& (z==149))||
	((y==149)&& (z==150))||
	((y==150)&& (z==149))||
	((y==150)&& (z==150))
	)
   {
		if (value != score -> end()) {
			//ofile2 << zz <<"  " <<*(value->second)/keV << G4endl;
			//ofile2 << xx << "  " << yy << "  " << zz <<"  " <<*(value->second)/keV << G4endl;
			if (xx<minX)minX=xx;
			if (xx>maxX)maxX=xx;
			G4double edep = *(value->second);///keV;
			voxelDose vD;
		    vD.xc=xx;vD.yc=yy;vD.zc=zz;//vD.edep=edep;
			vD.dose = edep;//edep/voxelMass;
		    //profileZmm.push_back(vD);
			profileX.push_back(vD);
//G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep = "<<G4BestUnit(vD.edep, "Energy")<<"; Dose = "
				//<<G4BestUnit(vD.dose,"Dose")<<G4endl;
		}
	}
	*/
//============================================================================
      if (value != score -> end()) 
        {
        // Print in the ASCII output file the information
        ofile << xx << "  " << yy << "  " << zz <<"  " <<*(value->second)/keV << G4endl;	
           
        // Save the same information in the output analysis file
        analysisManager = G4AnalysisManager::Instance();
        analysisManager -> FillNtupleDColumn(0, xx);
        analysisManager -> FillNtupleDColumn(1, yy);
        analysisManager -> FillNtupleDColumn(2, zz);
        analysisManager -> FillNtupleDColumn(3, *(value->second)/keV);
        analysisManager -> AddNtupleRow();
        }
     }//z 
    } //y
  } //x
  
 
  //==================================================
  int nmm;
 /* nmm = profileZmm.size();
  G4cout<<"Zmm size: "<<nmm<<"; minZ = "<<minZ<<"; maxZ = "<<maxZ<<G4endl;
  for(int i = 0; i < nmm; i++) {
	 voxelDose vD = profileZmm[i];

	 G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
  }
  nmm = profileZmp.size();
  G4cout<<"Zmp size: "<<nmm<<G4endl;
  for(int i = 0; i < nmm; i++) {
	 voxelDose vD = profileZmp[i];

	 G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
  }
  nmm = profileZpm.size();
  G4cout<<"Zpm size: "<<nmm<<G4endl;
  for(int i = 0; i < nmm; i++) {
	 voxelDose vD = profileZpm[i];

	 G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
  }
  nmm = profileZpp.size();
  G4cout<<"Zpp size: "<<nmm<<G4endl;
  for(int i = 0; i < nmm; i++) {
	 voxelDose vD = profileZpp[i];

	 G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
  }*/
 //nmm = profileZ.size();
  //G4cout<<"Z size: "<<nmm<<G4endl;
  //for(int i = 0; i < nmm; i++) {
	// voxelDose vD = profileZ[i];

	 //G4cout<<"; xc = "<<vD.xc<<"; yc = "<<vD.yc<<"; zc = "<<vD.zc<<"; Edep[keV] = "<<vD.edep<<G4endl;
  //}
  
/*  nmm= profileZ.size();
  G4double number = 4.0;//4 beans around center  
  G4double sum = 0.0;
  //double conversion = 1000.0*1000.0*NumberOfQuantaEmitted/NumberOfQuantaSimulated;//micro
  double conversion = NumberOfQuantaEmitted/NumberOfQuantaSimulated;
  //testing sum of squares:
  //conversion=conversion*conversion/gray;//WORKS!!!!!!!!!!!!!!
  if(iCall==1){
	  //print only once!
  G4cout<<"Radioactive source is placed in the center at (x,y,z) = (0,0,0) "<<G4endl;
  G4cout<<"Energy deposition (EDEP) is scored in each voxel (1mm size) of the mesh"<<G4endl;
  //G4cout<<"Voxel center coordinate on a particular central axis and EDEP profile (keV): "<<G4endl;
  G4cout<<"Voxel center coordinate on a particular central axis and DOSE profile (Gy): "<<G4endl;
  } else if (iCall==2){
	  conversion=conversion*conversion/gray;//never go into here
  }
  G4cout<<"z-axis"<<G4endl;
  G4cout<<"BEGIN_DATA"<<G4endl;
  for (G4double z = minZ; z<=maxZ; z=z+1.0){//1 mm increment
	sum = 0.0;
	for(int i = 0; i < nmm; i++) {
		voxelDose vD = profileZ[i];
		if (vD.zc==z) sum = sum + vD.dose;//vD.edep;
	}
	sum = sum /number;
	sum=sum*conversion/gray;//Gy->conversion
G4cout<<z<<" "<<sum<<G4endl;
	centerZ.push_back(z);
	finalProfileZ.push_back(sum);
  }
  G4cout<<"END_DATA"<<G4endl;

  nmm= profileY.size();
  sum = 0.0;
  G4cout<<"y-axis"<<G4endl;
  G4cout<<"BEGIN_DATA"<<G4endl;
  for (G4double y = minY; y<=maxY; y=y+1.0){//1 mm increment
	sum = 0.0;
	for(int i = 0; i < nmm; i++) {
		voxelDose vD = profileY[i];
		if (vD.yc==y) sum = sum + vD.dose;//vD.edep;
	}
	sum = sum /number;
	sum=sum*conversion/gray;//Gy->conversion
G4cout<<y<<" "<<sum<<G4endl;
	centerY.push_back(y);
	finalProfileY.push_back(sum);
  }
  G4cout<<"END_DATA"<<G4endl;

   nmm= profileX.size();
  sum = 0.0;
  G4cout<<"x-axis"<<G4endl;
  G4cout<<"BEGIN_DATA"<<G4endl;
  for (G4double x = minX; x<=maxX; x=x+1.0){//1 mm increment
	sum = 0.0;
	for(int i = 0; i < nmm; i++) {
		voxelDose vD = profileX[i];
		if (vD.xc==x) sum = sum + vD.dose;//vD.edep;
	}
	sum = sum /number;
	sum=sum*conversion/gray;//Gy->conversion
G4cout<<x<<" "<<sum<<G4endl;
	centerX.push_back(x);
	finalProfileX.push_back(sum);
  }
  G4cout<<"END_DATA"<<G4endl;*/
  //========================================================
  
  ofile << std::setprecision(6);//ofile2 << std::setprecision(6);

  // Close the output ASCII file
  ofile.close();//ofile2.close();

  // Close the output brachytherapy.root 
  analysisManager -> Write();
  analysisManager -> CloseFile();

  //clear for reuse
  ///centerZ.clear(); finalProfileZ.clear();centerY.clear(); finalProfileY.clear();centerX.clear(); finalProfileX.clear();
  //profileZ.clear();profileX.clear();profileY.clear();
}

