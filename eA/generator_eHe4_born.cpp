#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TVectorT.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unistd.h>

#include "ep_constants.h"

const double mHe4=3.7273794118; // Mass of He-4 nucleus in GeV

using namespace std;

double Gch(double QSq)
{
  // Parameterization from Frosch et al., Phys. Rev. 160, no. 4, p. 874 (1967)
  const double a=0.316/hbarc; // converted from fm to GeV^-1
  const double b=0.681/hbarc; // converted from fm to GeV^-1
  const double n=6.0;

  return (1.-pow(a*a*QSq,n))*exp(-b*b*QSq);
}

double dSigdOmega(double E1, double cosTheta)
{
  const double E3 = E1*mHe4/(mHe4+E1*(1.-cosTheta));
  const double cosSqThetaOver2=0.5*(1.+cosTheta);
  const double sinSqThetaOver2=0.5*(1.-cosTheta);

  const double QSq = 2.*mHe4*(E1-E3);
  const double sigma_Mott = nbGeVSq*alpha*alpha*E3*cosSqThetaOver2/(4.*E1*E1*E1*sinSqThetaOver2*sinSqThetaOver2);
  const double G=Gch(QSq);
  
  return sigma_Mott * G*G;
}

const int requiredArgs=4;

int main(int argc, char ** argv)
{
  if (argc < requiredArgs)
    {
      cerr << "Wrong number of arguments. Usage is:\n"
           << "    generator /path/to/output/file [Nevents] [Beam energy (GeV)] [optional arguments...]\n"
           << "        Optional arguments:\n"
           << "            -t: minimum theta [deg.]\n"
           << "            -T: maximum theta [deg.]\n"
           << "\n\n";

      return -1;
    }

  // Set up the output file and the tree
  const int Nevents=atoi(argv[2]);
  TFile * outfile = new TFile(argv[1],"RECREATE");
  TTree * outtree = new TTree("T","e-He4 generator tree");
  const double E1=atof(argv[3]);

  double cosTheta_e_max = TMath::Cos(1.*M_PI/180.);
  double cosTheta_e_min = TMath::Cos(179.*M_PI/180.);

  // optional arguments
  int c;
  while ((c = getopt (argc-requiredArgs+1, &argv[requiredArgs-1], "t:T:")) != -1)
    {
      switch (c)
        {
        case 't':
          cosTheta_e_max = cos(atof(optarg)*M_PI/180.);
          break;
        case 'T':
          cosTheta_e_min = cos(atof(optarg)*M_PI/180.);
          break;
	default:
          cerr << "Optional argument not found. Aborting...\n\n";
          return -1;
        }
    }

  // Now that the generator settings have been established, write them out to a TVectorT
  TVectorT<double> runInfo(4);
  runInfo[0]=Nevents;
  runInfo[1]=E1;
  runInfo[2]=cosTheta_e_max;
  runInfo[3]=cosTheta_e_min;
  runInfo.Write("runInfo");
  
  // Set up a random number generator
  TRandom3 * myRand = new TRandom3(0); 
 
  // Tree info
  double mom_e,theta_e,phi_e,mom_He, theta_He, phi_He, weight;
  outtree->Branch("mom_e",&mom_e, "mom_e/D");
  outtree->Branch("theta_e",&theta_e, "theta_e/D");
  outtree->Branch("phi_e",&phi_e, "phi_e/D");
  outtree->Branch("mom_He",&mom_He, "mom_He/D");
  outtree->Branch("theta_He",&theta_He, "theta_He/D");
  outtree->Branch("phi_He",&phi_He, "phi_He/D");
  outtree->Branch("weight",&weight, "weight/D");
  
  // Loop over events
  for (int event=0 ; event < Nevents ; event++)
    {
      if (event %100000==0)
	{
	  cerr << "Working on event " << event << " out of " << Nevents << "...\n";
	}
      
      // First step is to generate an electron scattering angle
      double cosTheta_e = cosTheta_e_min + (cosTheta_e_max - cosTheta_e_min)*myRand->Rndm();
      theta_e = TMath::ACos(cosTheta_e);
      
      // This can be used to establish kinematics
      mom_e=E1*mHe4 / (mHe4+E1*(1.-cosTheta_e));
      phi_e = myRand->Rndm()*2.*M_PI - M_PI;

      TVector3 p1(0.,0.,E1);
      TVector3 p3(mom_e*TMath::Sin(theta_e)*TMath::Cos(phi_e), mom_e*TMath::Sin(theta_e)*TMath::Sin(phi_e), mom_e*cosTheta_e);
      TVector3 p4 = p1 - p3;
      theta_He = p4.Theta();
      phi_He = p4.Phi();
      mom_He = p4.Mag();

      //Weight is: d^N sig / P(cosTheta, phi)
      weight = dSigdOmega(E1, cosTheta_e) * (2.*M_PI) *  (cosTheta_e_max - cosTheta_e_min);

      // Sanitize
      if (weight<0.)
	weight=0.;
      
      outtree->Fill();
    }

  outfile->Write();
  outfile->Close();

  return 0;
}
