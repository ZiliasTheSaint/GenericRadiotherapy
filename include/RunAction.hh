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

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include <vector> //just in case for below vectors
class G4Run;
class DetectorConstruction;

class RunAction : public G4UserRunAction
{
public:
  RunAction();
  ~RunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run* );
  
  void Fill(G4int index, G4double doseDeposit);
  void FillEnergy(G4int index, G4double energyDeposit);

  //for better statistics (covariance built in)
  void FillEnergyTotal(G4double eventEnergy);

  //variance x2/(n) - mean^2 can be negative due to roundoff errors!
  G4double basicAnalyze(G4double x, G4double x2, G4int n){
	  /*G4double temp=x;G4double temp2=x2;G4int nevents=n;
	  temp=temp/nevents;temp2=temp2/nevents;
	  temp2 = temp2 - temp*temp;
	  if (nevents>1) temp2=temp2/(nevents-1.0);
	  if (temp2 >0.)
	   temp2 = sqrt(temp2); 

	  //=========percent
	  if (temp!=0.0){
	   temp2 = std::min(100.0*temp2/temp,99.9);
      } 
	  return temp2;*/
	  if (n < 1 || x == 0.0) return 0.0;
	  
	  const G4double mean = x/n;
      G4double var = x2/n - mean*mean;
      
	  if (n > 1) var /= (n - 1.0);
      
	  const G4double relErr = (var > 0.0) ? 100.0*std::sqrt(var)/std::fabs(mean) : 0.0;//in %!
      return std::min(relErr, 99.9);//in %!!
  }
private:
	void totalRunEnergyDeposit(G4int n);
	DetectorConstruction*    phantom;     //pointer to the geometry
	std::map<G4int,G4double> dose;
	std::map<G4int,G4double> dose2;
	std::map<G4int,G4double> EDEP;
	std::map<G4int,G4double> EDEP2;

	G4double totalEdepRun;
	G4double totalEdepRun2;

	double minRelErrX,maxRelErrX,minRelErrY,maxRelErrY,minRelErrZ,maxRelErrZ;

	//std::vector<double> projyProfile;std::vector<double> ycoordProfile;
	//std::vector<double> projxProfile;std::vector<double> xcoordProfile;
	//std::vector<double> projzProfile;std::vector<double> zcoordProfile;

	std::vector<G4double> projyProfile;std::vector<G4double> ycoordProfile;
	std::vector<G4double> projxProfile;std::vector<G4double> xcoordProfile;
	std::vector<G4double> projzProfile;std::vector<G4double> zcoordProfile;

	void PrintAxisProfile(const char* axisName, G4int nBins, G4double voxelWidth,
                        const std::vector<G4double>& sum, const std::vector<G4double>& sum2,
                        G4int nVoxelsInBlock, G4int nEvents, G4double localConversion
						//added						
						,double& min, double& max
						//not really needed
						,std::vector<G4double>& coordOut, std::vector<G4double>& doseOut
						);
};
#endif



