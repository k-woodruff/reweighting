#include <cmath>
#include <iostream>

#include "ncformfactor.h"

double P_EFF(double q2);
double N_EFF(double q2);
double P_MFF(double q2);
double N_MFF(double q2);

double pi  =  3.1415927;
/* neutron magnetic moment--units of  */
double muN = -1.91304;
/* vector mass--units of GeV/c^2 */
double mv  = 0.840;
/* proton magnetic moment--units of  */
double muP = 2.7930;

/*Mass of proton in GeV*/
double MP  =  0.9382721;
/*Mass of muon in GeV*/
double Mmu =  0.105658;

/* (G_FERMI * hbar * c)^2, in units of 10^{-12} (fm/GeV)^2 */
double gfhbc2        =  5.2971;
/*Axial coupling constant*/
double g_a           =  1.2670;

//  This character allows you to select the parametrization you want
//  There is NO DEFAULT -- you must set a valid value for FF
//  G = GKex                 D = Dipole-Galster      K = Kelly
//  F = Friedrich-Walcher    A = Arrington-Sick      T = A-S with TPE correction
//  B = BBA05 (matches genie CCQE)
char FF  = 'B';

/*anomalous magnetic moments*/ 
double k_rho = 5.51564; 
double k_omega = 0.4027;
double krho_prime = 12.0; 
double komega_prime = -2.973; 
double k_phi = 0.01; 
double k_s = -0.120195; // = k_p + k_n
double kv = 3.705889; // = k_p - k_n
double mu_phi = 0.2; 
double alpha1 = 0.0781808;
double alpha2 = 0.0632907; 
/*fixed masses in GeV*/
double lambda_d = 1.181; 
double lambda_qcd = 0.150;
double lambda1 = 0.93088; 
double lambda2 = 2.6115;  
double m_rho = 0.776; 
double m_omega = 0.784; 
double mrho_prime = 1.45; 
double del1 = 0.03465; 
double del2 = 0.04374;
double m1_rho = 0.741; //  = m_rho - del1
double m2_rho = 0.732; // = m_rho - del2
double momega_prime = 1.419; 
double m_phi = 1.019; 
double gf_rho = 0.5596; 
double gf_omega = 0.7021; 
double gfprime_rho = 0.0072089; 
double gfprime_omega = 0.164; 
double gf_phi = -0.1711;
/*momentum cutoff in (GeV / C)^2*/    
double Q2_1 = 0.3176; 
double Q2_2 = 0.1422; 

//  Continued-fraction function used in Arrington-Sick model
double GCF(double q2, double b1, double b2, double b3, double b4, double b5)
{ return 1/(1+b1*q2/(1+b2*q2/(1+b3*q2/(1+b4*q2/(1+b5*q2))))); }


// These functions calculate various form factor parametrizations
double Qt(double q2)
{ return q2 * log((pow(lambda_d,2) + q2) / pow(lambda_qcd,2)) / log(pow(lambda_d,2) / pow(lambda_qcd,2)); } 
double fem(double mx, double q2)
{ return (pow(mx,2) / (pow(mx,2) + q2)); }
double f(double q2, double lambda)
{ return pow(lambda,2) / (pow(lambda,2) + q2); }
double f1_had(double q2)
{ return f(Qt(q2),lambda1) * f(Qt(q2),lambda2); } 
double f2_had(double q2)
{ return f(Qt(q2),lambda1) * pow(f(Qt(q2),lambda2),2); } 
double fs1_had(double q2)
{ return f1_had(q2) * pow((q2 / (pow(lambda1,2) + q2)),1.5); } 
double fs2_had(double q2)
{ return f2_had(q2) * pow(((pow(mu_phi,2) + q2) / pow(mu_phi,2)) * (pow(lambda1,2) / (pow(lambda1,2) + q2)),1.5); } 
double f1_pQCD(double q2)
{ return (f(Qt(q2),lambda_d)) * (f(Qt(q2),lambda2)); } 
double f2_pQCD(double q2)
{ return f(Qt(q2),lambda_d) * pow(f(Qt(q2),lambda2),2);  } 

/*  Definition of GKex model via defined quantities*/

double F0_1(double q2)
{ return gf_omega * fem(m_omega, q2) * f1_had(q2) + 
  gfprime_omega * fem(momega_prime,q2)*f1_had(q2) + 
  gf_phi * fem(m_phi,q2) * fs1_had(q2) + 
  (1 - gf_omega - gfprime_omega) * f1_pQCD(q2); } 

double F0_2(double q2)
{ return k_omega * gf_omega * fem(m_omega,q2)*f2_had(q2) + 
  komega_prime * gfprime_omega * fem(momega_prime,q2)*f2_had(q2)+
  k_phi * gf_phi * fem(m_phi, q2) * fs2_had(q2) + 
  (k_s - (k_omega * gf_omega) - (komega_prime * gfprime_omega) - (k_phi *gf_phi)) * f2_pQCD(q2); } 

double F1_1(double q2)
{ return gf_rho * fem(m1_rho, q2) * f1_had(q2)  * 
  ((1 - alpha1) + (alpha1 / pow(1 + q2 / Q2_1,2))) + 
  gfprime_rho * fem(mrho_prime,q2) * f1_had(q2) +
  (1 - gf_rho - gfprime_rho) * f1_pQCD(q2) ; } 

double F1_2(double q2)
{ return k_rho * gf_rho * fem(m2_rho,q2) * f2_had(q2) * 
  ((1 - alpha2) + (alpha2 / (1 + q2 / Q2_2))) + 
  krho_prime * gfprime_rho * fem(mrho_prime,q2) * f2_had(q2) + 
  (kv - (k_rho*gf_rho) - (krho_prime * gfprime_rho)) * f2_pQCD(q2); }
    
//cross section for charged current
double sigplus_CC(double Ev, double q2, double MA,
                  int nusign)
{

  //Cos^{2}(theta_cabbibo)                                          
  double cos2_cabbibo  =  0.94711824;  
  double tau = q2 / pow(2.*MP,2);

  double FP1  = (P_EFF(q2) + tau*P_MFF(q2)) / (1+tau);
  double FP2  = (P_MFF(q2) - P_EFF(q2)) / (1+tau);
  double FN1  = (N_EFF(q2) + tau*N_MFF(q2))/(1+tau); 
  double FN2  = (N_MFF(q2) - N_EFF(q2)) / (1+tau);

  double GVA  = (-1.)*(g_a/pow(1+q2/pow(MA,2),2));
  double FV1  = (FP1 - FN1);
  double FV2  = (FP2 - FN2);

  double W_CC = 4*((Ev/MP) - tau - (pow(Mmu,2)/(4*pow(MP,2))));
  double A_CC = ((pow(Mmu,2) + q2) / 4)*((pow(GVA,2)) * (1+tau) -(pow(FV1,2) - (tau*pow(FV2,2)))*(1-tau) + (4*tau*FV1*FV2));
  double B_CC = -(q2/4)*GVA*(FV1 + FV2);
  double C_CC = (q2/(64*tau)) * (pow(GVA,2) + pow(FV1,2) + (tau*pow(FV2,2)));

  double GeVfm = 1.e-12/pow(0.19732697,2); //1/hbarc^2
	return (gfhbc2*cos2_cabbibo/(2*pi))*(1/pow(Ev,2))*(A_CC + (nusign)*B_CC*W_CC + C_CC*pow(W_CC,2)) * GeVfm;

} 
	

//cross section for neutral current
double sigplus(double Ev, double q2, char gas, 
               double delta_s, double MA,
               double lambda_a, double s_a,
               double mu_s, double rho_s,
               int nusign, int nucsign)
{ 
  /*
  // fitting parameters for strange form factors
  // default values based on Summer 2016 fit, Modified Dipole Model
  rho_s    = -0.10; 
  mu_s     =  0.056;
  delta_s  = -0.29;
  lambda_a =  1.1;
  s_a      =  0.4;
  // Axial mass constant
  MA       =  1.001;
  */

  /* sin^2(Theta_W), weak mixing angle */
  double sin2_weinberg =  0.2377;
  double tau = q2 / pow(2.*MP,2);

  double GAS; 
  switch(gas){ 
    case 'a' :
      GAS = (delta_s + s_a*q2)/pow((1+q2/pow(lambda_a, 2)),2);
      break;   
  
    case 'b' :
      GAS = delta_s /( 1 + (lambda_a*q2)/(1+(s_a*q2)) ); 
      break;
  
    case 'c' :
      GAS = delta_s /pow(1+q2/pow(lambda_a,2),2); 
      break;
  
    default :
        throw "::Invalid form factor (gas) model.";
  }

  double FP1 = (P_EFF(q2) + tau*P_MFF(q2)) / (1+tau);
  double FP2 = (P_MFF(q2) - P_EFF(q2)) / (1+tau);
  double FN1 = (N_EFF(q2) + tau*N_MFF(q2))/(1+tau) ; 
  double FN2 = (N_MFF(q2) - N_EFF(q2)) / (1+tau);

  double GACC = g_a/pow(1+q2/pow(MA,2),2);
  double GZA  = (nucsign)*(-0.5)*GACC + 0.5*GAS;

  double FS1 = (rho_s*tau + tau*mu_s) / (1+ tau);
  double FS2 = (-rho_s*tau + mu_s) / (1+tau); 
  //double FZ1 = 0.5*((1-4*sin2_weinberg)*(FP1 - FN1 - FS1));
  double FZ1,FZ2;
  if(nucsign == -1)
  {
    FZ1 = 0.5*((1-4*sin2_weinberg)*(FN1 - FP1 - FS1));
    FZ2 = 0.5*((1-4*sin2_weinberg)*FN2 - FP2 - FS2);
  }
  else
  {
    FZ1 = 0.5*((1-4*sin2_weinberg)*(FP1 - FN1 - FS1));
    FZ2 = 0.5*((1-4*sin2_weinberg)*FP2 - FN2 - FS2);
  }

  double W = 4.*(Ev/MP - tau);
  double A = (q2/4.) * ((pow(GZA,2)*(1+tau)) - (pow(FZ1,2) - tau*pow(FZ2,2))*(1-tau) + 4*tau*FZ1*FZ2);
  double B = -(q2/4.*GZA)*(FZ1 + FZ2);
  double C = q2/(64.*tau) * (pow(GZA,2)+pow(FZ1,2) + tau*(pow(FZ2,2))); 

  double GeVfm = 1.e-12/pow(0.19732697,2); //1/hbarc^2
  return (gfhbc2/(2*pi))*(1/pow(Ev,2))*(A + (nusign)*B*W + C*pow(W,2) ) * GeVfm; 
} 

double P_EFF(double q2) {
  double denom;
  
  double tau = q2 / pow(2.*MP,2);
  double a0=1.0;
  double a1=-0.24;
  double b1=10.98;
  double b2=12.82;
  double b3=21.97;
  
  double a10 = 1.041;
  double a11 = 0.765;
  double a20 = -0.041;
  double a21 = 6.2;
  double ab = -0.23;
  double qb = 0.07;
  double sigb = 0.27;    
  double q;

  // parameters for A-S model 
  double b1_gep = 3.440;
  double b2_gep = -0.178;
  double b3_gep = -1.212;
  double b4_gep = 1.176;
  double b5_gep = -0.284;
  
  double b1_gep_tpe = 3.478;
  double b2_gep_tpe = -0.140;
  double b3_gep_tpe = -1.311;
  double b4_gep_tpe = 1.128; 
  double b5_gep_tpe = -0.233;

  // parameters for BBA05 model
  double a0_bb = 1.;
  double a1_bb = -0.0578;
  double a2_bb = 0.;
  double b1_bb = 11.1;
  double b2_bb = 13.6;
  double b3_bb = 33.0;
  double b4_bb = 0.;

  switch (FF) {

  case 'A' :
    return GCF(q2,b1_gep,b2_gep,b3_gep,b4_gep,b5_gep);
    break;

  case 'T' :
    return GCF(q2,b1_gep_tpe,b2_gep_tpe,b3_gep_tpe,b4_gep_tpe,b5_gep_tpe);
    break;


  case 'G' : 
    return (F0_1(q2) + F1_1(q2)) / 2 - tau*(F0_2(q2) + F1_2(q2)) / 2;
    break; 
    
  case 'D':
    denom = 1+q2/mv/mv;
    return 1/denom/denom;
    break;
    
  case 'K':
    return (a0 + a1*tau)
      /(1 + b1*tau + b2*tau*tau + b3*tau*tau*tau);
    break;
    
  case 'B':
    return (a0_bb + a1_bb*tau + a2_bb*tau*tau)
      /(1 + b1_bb*tau + b2_bb*tau*tau + b3_bb*tau*tau*tau + b4_bb*tau*tau*tau*tau);
    break;
    
  case 'F':
    q=sqrt(q2);  
    return a10/pow(1+q2/a11,2) + a20/pow(1+q2/a21,2)
      + ab*q2*( exp(-pow((q-qb)/sigb,2)/2) + exp(-pow((q+qb)/sigb,2)/2) );
    break;
 
  default :
      throw "::Invalid form factor (FF) model.";
  }
}

double N_EFF(double q2) {
  double denom;
  
  double A=1.70;
  double B=3.30;
  double lambda2=0.71;
  
  double tau = q2 / pow(2.*MP,2);
  double a10 = 1.04;
  double a11 = 1.73;
  double a20 = -1.04;
  double a21 = 1.54;
  double ab = 0.23;
  double qb = 0.29;
  double sigb = 0.20;
  double q;
  
  // parameters for A-S model 
  double b1_gen = 0.977;
  double b2_gen = -20.82;
  double b3_gen = 22.02;

  // parameters for BBA05 model
  double a0_bb = 0.;
  double a1_bb = 1.25;
  double a2_bb = 1.3;
  double b1_bb = -9.86;
  double b2_bb = 305.;
  double b3_bb = -758.;
  double b4_bb = 802.;
  
  switch (FF) {

  case 'A' :
  case 'T' :
    return 0.484*q2*GCF(q2,b1_gen,b2_gen,b3_gen,0.0,0.0);
    break;

  case 'G' :
    return (F0_1(q2) - F1_1(q2)) / 2 - tau*(F0_2(q2) - F1_2(q2)) / 2; 
    break; 
    
  case 'D':
    denom = 1+q2/mv/mv;
    return -tau*muN/denom/denom/(1+5.6*tau);
    break;
    
  case 'K':
    return (A*tau)/(1 + B*tau)/(1 + q2/lambda2)/(1 + q2/lambda2);
    break;
    
  case 'B':
    return (a0_bb + a1_bb*tau + a2_bb*tau*tau)
      /(1 + b1_bb*tau + b2_bb*tau*tau + b3_bb*tau*tau*tau + b4_bb*tau*tau*tau*tau);
    break;
    
  case 'F':
    q=sqrt(q2);
    return a10/pow(1+q2/a11,2) + a20/pow(1+q2/a21,2)
      + ab*q2*( exp(-pow((q-qb)/sigb,2)/2) + exp(-pow((q+qb)/sigb,2)/2) );
    break;

  default :
      throw "::Invalid form factor (FF) model.";
  }  
}


double P_MFF(double q2) {
  double denom;
  
  double tau = q2 / pow(2.*MP,2);
  double a0=1.0;
  double a1=0.12;
  double b1=10.97;
  double b2=18.86;
  double b3=6.55;
  
  double a10 = 1.002;
  double a11 = 0.749;
  double a20 = -0.002;
  double a21 = 6.0;
  double ab = -0.13;
  double qb = 0.35;
  double sigb = 0.21;
  double q;
  
  // parameters for A-S model 
  double b1_gmp = 3.173;
  double b2_gmp = -0.314;
  double b3_gmp = -1.165;
  double b4_gmp = 5.619;
  double b5_gmp = -1.087;
  
  double b1_gmp_tpe = 3.224;
  double b2_gmp_tpe = -0.313;
  double b3_gmp_tpe = -0.868;
  double b4_gmp_tpe = 4.278;
  double b5_gmp_tpe = -1.102;

  // parameters for BBA05 model
  double a0_bb = 1.;
  double a1_bb = 0.15;
  double a2_bb = 0.;
  double b1_bb = 11.1;
  double b2_bb = 19.6;
  double b3_bb = 7.54;
  double b4_bb = 0.;

  switch (FF) {

  case 'A' :
    return muP*GCF(q2,b1_gmp,b2_gmp,b3_gmp,b4_gmp,b5_gmp);
    break;

  case 'T' :
    return muP*GCF(q2,b1_gmp_tpe,b2_gmp_tpe,b3_gmp_tpe,b4_gmp_tpe,b5_gmp_tpe);
    break;

  case 'G': 
    return (F0_1(q2) + F1_1(q2)) / 2 + (F0_2(q2) + F1_2(q2)) / 2; 
    break; 
    
  case 'D':
    denom = 1+q2/mv/mv;
    return muP/denom/denom;
    break;
    
  case 'K':
    return muP*(a0 + a1*tau)
      /(1 + b1*tau + b2*tau*tau + b3*tau*tau*tau);
    break;
    
  case 'B':
    return muP*(a0_bb + a1_bb*tau + a2_bb*tau*tau)
      /(1 + b1_bb*tau + b2_bb*tau*tau + b3_bb*tau*tau*tau + b4_bb*tau*tau*tau*tau);
    break;
    
  case 'F':
    q=sqrt(q2);
    return muP*a10/pow(1+q2/a11,2) + a20/pow(1+q2/a21,2)
      + ab*q2*( exp(-pow((q-qb)/sigb,2)/2) + exp(-pow((q+qb)/sigb,2)/2) );
    break;

  default :
    throw "::Invalid form factor (FF) model.";
  }
}

double N_MFF(double q2) {
 
  double denom;

  double tau = q2 / pow(2.*MP,2);
  double a0=1.0;
  double a1=2.33;
  double b1=14.72;
  double b2=24.20;
  double b3=84.1;

  double a10 = 1.012;
  double a11 = 0.770;
  double a20 = -0.012;
  double a21 = 6.8;
  double ab = -0.28;
  double qb = 0.33;
  double sigb = 0.14;
  double q;

  // parameters for A-S model 
  double b1_gmn = 3.297;
  double b2_gmn = -0.258;
  double b3_gmn = 0.001;
  
  // parameters for BBA05 model
  double a0_bb = 1.;
  double a1_bb = 1.81;
  double a2_bb = 0.;
  double b1_bb = 14.1;
  double b2_bb = 20.7;
  double b3_bb = 68.7;
  double b4_bb = 0.;
  
  switch (FF) {

  case 'A' :
  case 'T' :
    return muN*GCF(q2,b1_gmn,b2_gmn,b3_gmn,0.0,0.0);
    break;

  case 'G': 
    return (F0_1(q2) - F1_1(q2)) / 2 + (F0_2(q2) - F1_2(q2)) / 2; 
    break; 
    
  case 'D': 
    denom = 1+q2/mv/mv;
    return muN/denom/denom;
    break;
    
  case 'K':
    return muN*(a0 + a1*tau)
      /(1 + b1*tau + b2*tau*tau + b3*tau*tau*tau);
    break;
    
  case 'B':
    return muN*(a0_bb + a1_bb*tau + a2_bb*tau*tau)
      /(1 + b1_bb*tau + b2_bb*tau*tau + b3_bb*tau*tau*tau + b4_bb*tau*tau*tau*tau);
    break;
    
  case 'F':
    q=sqrt(q2);
    return muN*a10/pow(1+q2/a11,2) + a20/pow(1+q2/a21,2)
      + ab*q2*( exp(-pow((q-qb)/sigb,2)/2) + exp(-pow((q+qb)/sigb,2)/2) );
    break;

  default :
    throw "::Invalid form factor (FF) model.";
  }
}

