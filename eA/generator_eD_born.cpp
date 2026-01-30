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

using namespace std;

const double mD=1.875612945;  // Mass of the deuteron (deuterium nucleus) in GeV

double GAbbott1(double QSq, double G0, double Q0, double a1, double a2, double a3, double a4, double a5)
{
  const double QSqInvFM = QSq/(hbarc*hbarc);
  // This parameterization comes from D. Abbott et al., Eur. Phys. J. A 7, 421 (2000)
  return G0 * (1. - QSqInvFM/(Q0*Q0)) * (1. + a1*QSqInvFM + a2*QSqInvFM*QSqInvFM + a3*QSqInvFM*QSqInvFM*QSqInvFM +
  					 a4*QSqInvFM*QSqInvFM*QSqInvFM*QSqInvFM + a5*QSqInvFM*QSqInvFM*QSqInvFM*QSqInvFM*QSqInvFM);
}

// Parameter values come from J. Zhou et al., Eur. Phys. J. A (2023) 59:256
// https://doi.org/10.1140/epja/s10050-023-01174-6
// Parameter values are in units of powers of fm^-1
double GC(double QSq)
{
  return GAbbott1(QSq,1.,4.21, 0.674, 0.02246,9.806e-3,-2.709E-4,3.793E-6);
}

double GM(double QSq)
{
  return GAbbott1(QSq,1.714,7.37,0.5804,0.08701,-3.624E-3,3.448E-4,-2.818E-6);
}

double GQ(double QSq)
{
  return GAbbott1(QSq,25.83,8.10,0.876,-5.656e-2,1.933e-2,-6.734E-4,9.348e-6);
}

double dSigdOmega(double E1, double cosTheta)
{
  const double E3 = E1*mD/(mD+E1*(1.-cosTheta));
  const double cosSqThetaOver2=0.5*(1.+cosTheta);
  const double sinSqThetaOver2=0.5*(1.-cosTheta);
  const double tanSqThetaOver2=sinSqThetaOver2/cosSqThetaOver2;

  const double QSq = 2.*mD*(E1-E3);
  const double eta = QSq/(4.*mD*mD);

  const double A=GC(QSq)*GC(QSq) + 8./9.*eta*eta*GQ(QSq)*GQ(QSq) + 2./3.*eta*GM(QSq)*GM(QSq);
  const double B=4./3.*eta*(1.+eta)*GM(QSq)*GM(QSq);

  const double sigma_Mott = nbGeVSq*alpha*alpha*E3*cosSqThetaOver2/(4.*E1*E1*E1*sinSqThetaOver2*sinSqThetaOver2);
  return sigma_Mott * (A + B*tanSqThetaOver2);
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
  double mom_e,theta_e,phi_e,mom_D, theta_D, phi_D, weight;
  outtree->Branch("mom_e",&mom_e, "mom_e/D");
  outtree->Branch("theta_e",&theta_e, "theta_e/D");
  outtree->Branch("phi_e",&phi_e, "phi_e/D");
  outtree->Branch("mom_D",&mom_D, "mom_D/D");
  outtree->Branch("theta_D",&theta_D, "theta_D/D");
  outtree->Branch("phi_D",&phi_D, "phi_D/D");
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
      mom_e=E1*mD / (mD+E1*(1.-cosTheta_e));
      phi_e = myRand->Rndm()*2.*M_PI - M_PI;

      TVector3 p1(0.,0.,E1);
      TVector3 p3(mom_e*TMath::Sin(theta_e)*TMath::Cos(phi_e), mom_e*TMath::Sin(theta_e)*TMath::Sin(phi_e), mom_e*cosTheta_e);
      TVector3 p4 = p1 - p3;
      theta_D = p4.Theta();
      phi_D = p4.Phi();
      mom_D = p4.Mag();

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
