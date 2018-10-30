#include <cmath>
#include <iostream>

#include "zeformfactor.h"

double tau;
double pi  =  3.1415927;
/* neutron magnetic moment--units of  */
double muN = -1.91304;
/* vector mass--units of GeV/c^2 */
double mv  = 0.840;
/* proton magnetic moment--units of  */
double muP = 2.7930;

/*Mass of proton in GeV*/
double MP  =  0.9382721;
/*Mass of neutron in GeV*/
double MN  =  0.9395656;
/*Mass of muon in GeV*/
double Mmu =  0.105658;
/*Mass of charged pion in GeV*/
double Mpi =  0.13957;

/* (G_FERMI * hbar * c)^2, in units of 10^{-12} (fm/GeV)^2 */
double gfhbc2        =  5.2971;
/*Axial coupling constant*/
double g_a           =  1.2670;

/* fixed parameters for z-expansion */
double tcut    = 4.*pow(Mpi,2);
double t0      = -0.7;
double tcut_ga = 9.*pow(Mpi,2);

//cross section for charged current
double sigplus_CC(double Ev, double q2, int nusign )
{

  //Cos^{2}(theta_cabbibo)                                          
  double cos2_cabbibo  =  0.94711824;  
  double tau = q2 / pow(2.*MP,2);

  double GVE = GEp(q2) - GEn(q2);
  double GVM = GMp(q2) - GMn(q2);

  double GA = GACC(q2);
  double FV1 = GVE;
  double FV2 = GVM;

  double W_CC = 4*((Ev/MP) - tau - (pow(Mmu,2)/(4*pow(MP,2))));
  double A_CC = ((pow(Mmu,2) + q2) / 4)*((pow(GA,2)) * (1+tau) -(pow(FV1,2) - (tau*pow(FV2,2)))*(1-tau) + (4*tau*FV1*FV2));
  double B_CC = (q2/4)*GA*(FV1 + FV2);
  double C_CC = (q2/(64*tau)) * (pow(GA,2) + pow(FV1,2) + (tau*pow(FV2,2)));

  double GeVfm = 1.e-12/pow(0.19732697,2); //1/hbarc^2
	return (gfhbc2*cos2_cabbibo/(2*pi))*(1/pow(Ev,2))*(A_CC - (nusign)*B_CC*W_CC + C_CC*pow(W_CC,2)) * GeVfm;

} 
	

//cross section for neutral current
double sigplus( double Ev, double q2,
                double delta_s, double MA_s, double a2_s,
                int nusign, int nucsign )
{ 

  /* sin^2(Theta_W), weak mixing angle */
  double sin2_weinberg =  0.2377;
  double tau = q2 / pow(2.*MP,2);

  double GVE, GVM;
  if(nucsign < 0)
  {
    GVE = -0.5*(GEp(q2) - GEn(q2)) - 2.*sin2_weinberg*GEn(q2);
    GVM = -0.5*(GMp(q2) - GMn(q2)) - 2.*sin2_weinberg*GMn(q2);
  }
  else
  {
    GVE = 0.5*(GEp(q2) - GEn(q2)) - 2.*sin2_weinberg*GEp(q2);
    GVM = 0.5*(GMp(q2) - GMn(q2)) - 2.*sin2_weinberg*GMp(q2);
  }

  double GZA = (nucsign)*0.5*GACC(q2) + 0.5*GAS(q2,delta_s,MA_s,a2_s);
  double FV1 = GVE;
  double FV2 = GVM;

  double W = 4.*(Ev/MP - tau);
  double A = (q2/4.) * ((pow(GZA,2)*(1+tau)) - (pow(FV1,2) - tau*pow(FV2,2))*(1-tau) + 4*tau*FV1*FV2);
  double B = (q2/4.)*GZA*(FV1 + FV2);
  double C = q2/(64.*tau) * (pow(GZA,2)+pow(FV1,2) + tau*(pow(FV2,2))); 

  double GeVfm = 1.e-12/pow(0.19732697,2); //1/hbarc^2
  return (gfhbc2/(2*pi))*(1/pow(Ev,2))*(A - (nusign)*B*W + C*pow(W,2) ) * GeVfm; 
} 

double GEp( double q2 )
{
  double z = (sqrt(tcut + q2) - sqrt(tcut - t0))/(sqrt(tcut + q2) + sqrt(tcut - t0));

  double a_ep[13];
  a_ep[0]  =  0.239163298067;
  a_ep[1]  = -1.109858574410;
  a_ep[2]  =  1.444380813060;
  a_ep[3]  =  0.479569465603;
  a_ep[4]  = -2.286894741870;
  a_ep[5]  =  1.126632984980;
  a_ep[6]  =  1.250619843540;
  a_ep[7]  = -3.631020471590;
  a_ep[8]  =  4.082217023790;
  a_ep[9]  =  0.504097346499;
  a_ep[10] = -5.085120460510;
  a_ep[11] =  3.967742543950;
  a_ep[12] = -0.981529071103;

  double gep = 0.;
  size_t alen = sizeof(a_ep)/sizeof(a_ep[0]);
  for(size_t k=0; k < alen; k++)
  {
    gep += a_ep[k]*pow(z,k);
  }

  return gep;
}

double GMp( double q2 )
{
  double z = (sqrt(tcut + q2) - sqrt(tcut - t0))/(sqrt(tcut + q2) + sqrt(tcut - t0));

  double a_mp[13];
  a_mp[0]  =  0.264142994136;
  a_mp[1]  = -1.095306122120;
  a_mp[2]  =  1.218553781780;
  a_mp[3]  =  0.661136493537;
  a_mp[4]  = -1.405678925030;
  a_mp[5]  = -1.356418438880;
  a_mp[6]  =  1.447029155340;
  a_mp[7]  =  4.235669735900;
  a_mp[8]  = -5.334045653410;
  a_mp[9]  = -2.916300520960;
  a_mp[10] =  8.707403067570;
  a_mp[11] = -5.706999943750;
  a_mp[12] =  1.280814375890;

  double gmp = 0.;
  size_t alen = sizeof(a_mp)/sizeof(a_mp[0]);
  for(size_t k=0; k < alen; k++)
  {
    gmp += a_mp[k]*pow(z,k);
  }

  return gmp*muP;
}

double GEn( double q2 )
{
  double z = (sqrt(tcut + q2) - sqrt(tcut - t0))/(sqrt(tcut + q2) + sqrt(tcut - t0));

  double a_en[11];
  a_en[0]  =  0.048919981379;
  a_en[1]  = -0.064525053912;
  a_en[2]  = -0.240825897382;
  a_en[3]  =  0.392108744873;
  a_en[4]  =  0.300445258602;
  a_en[5]  = -0.661888687179;
  a_en[6]  = -0.175639769687;
  a_en[7]  =  0.624691724461;
  a_en[8]  = -0.077684299367;
  a_en[9]  = -0.236003975259;
  a_en[10] =  0.090401973470;

  double gen = 0.;
  size_t alen = sizeof(a_en)/sizeof(a_en[0]);
  for(size_t k=0; k < alen; k++)
  {
    gen += a_en[k]*pow(z,k);
  }

  return gen;
}

double GMn( double q2 )
{
  double z = (sqrt(tcut + q2) - sqrt(tcut - t0))/(sqrt(tcut + q2) + sqrt(tcut - t0));

  double a_mn[11];
  a_mn[0]  =  0.257758326959;
  a_mn[1]  = -1.079540642058;
  a_mn[2]  =  1.182183812195;
  a_mn[3]  =  0.711015085833;
  a_mn[4]  = -1.348080936796;
  a_mn[5]  = -1.662444025208;
  a_mn[6]  =  2.624354426029;
  a_mn[7]  =  1.751234494568;
  a_mn[8]  = -4.922300878888;
  a_mn[9]  =  3.197892727312;
  a_mn[10] = -0.712072389946;

  double gmn = 0.;
  size_t alen = sizeof(a_mn)/sizeof(a_mn[0]);
  for(size_t k=0; k < alen; k++)
  {
    gmn += a_mn[k]*pow(z,k);
  }

  return gmn*muN;
}

double GACC( double q2 )
{
  double t0_ga = -0.28;
  double z = (sqrt(tcut_ga + q2) - sqrt(tcut_ga - t0_ga))/(sqrt(tcut_ga + q2) + sqrt(tcut_ga - t0_ga));

  double a_cc[9] = {-0.759,2.30,-0.6,-3.8,2.3,2.16,-0.896,-1.58,0.823};

  double gacc = 0.;
  size_t alen = sizeof(a_cc)/sizeof(a_cc[0]);
  for(size_t k=0;k < alen; k++)
  {
    gacc += a_cc[k]*pow(z,k);
  }

  return gacc;
}

double GAS( double q2, double a0, double a1, double a2 )
{
  double t0_s = -0.28;
  double z = (sqrt(tcut_ga + q2) - sqrt(tcut_ga - t0_s))/(sqrt(tcut_ga + q2) + sqrt(tcut_ga - t0_s));

  double a0_s = a0;
  double a1_s = a1;
  double a2_s = a2;

  double a3_s = -20.*a0_s - 10.*a1_s - 4.*a2_s;
  double a4_s =  45.*a0_s + 20.*a1_s + 6.*a2_s;
  double a5_s = -36.*a0_s - 15.*a1_s - 4.*a2_s;
  double a6_s =  10.*a0_s +  4.*a1_s + 1.*a2_s;

  return a0_s + a1_s*z + a2_s*z*z + a3_s*z*z*z 
       + a4_s*z*z*z*z + a5_s*z*z*z*z*z + a6_s*z*z*z*z*z*z;
}

/*
double GAS( double q2, double ds, double mas, double a2s )
{
  double t0_s = 0.;
  double z = (sqrt(tcut_ga + q2) - sqrt(tcut_ga - t0_s))/(sqrt(tcut_ga + q2) + sqrt(tcut_ga - t0_s));

  double a0_s = ds;
  double a1_s = -2.*ds/(mas*mas);
  double a2_s = a2s;

  return a0_s + a1_s*z + a2_s*z*z;
}
*/
