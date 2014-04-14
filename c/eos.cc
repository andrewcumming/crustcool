#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>
#include <stdlib.h>
#include "../h/root.h"
#include "../h/nr.h"

#define me 510.999
#define RADa 7.5657e-15

#include "../h/odeint.h"
#include "../h/eos.h"
#include <gsl/gsl_sf.h>

void* pt2Object;

// Wrappers for Potekhin's conductivity and EOS routines
extern "C"{
  void condegin_(double *temp,double *densi,double *B,double *Zion,double *CMI,
					double *CMI1,double *Zimp, double *RSIGMA,double *RTSIGMA,
					double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA);
  void eosmag_(double *Zion,double *CMI,double *RHO,double *TEMP,double *GAMAG,
     			double *DENS,double *GAMI,double *CHI,double *TPT,double *LIQSOL,
				double *PnkT,double *UNkT,double *SNk,double *CVE,double *CVI,double *CHIR,double *CHIT);
};


// ------------------------ initialise ----------------------------------

Eos::~Eos()
{
    free_vector(this->A,1,this->ns);
    free_vector(this->Z,1,this->ns);
    free_vector(this->X,1,this->ns);
}

Eos::Eos(int n)
{
	this->ns=n;
	this->A=vector(1,this->ns);
	this->Z=vector(1,this->ns);
	this->X=vector(1,this->ns);
	this->set_YZ2=0.0;
	this->Yn=0.0;
	this->set_Ye=0.0;
	this->set_Yi=0.0;
	this->accr=1;  // default crust composition is HZ with Fe56
	this->gamma_melt=175.0;
	this->Qimp=900.0; // treat the crust as a liquid for conductivities
	this->B=0.0; // default is unmagnetized 
	this->use_potek_cond = 1;
	this->use_potek_eos = 0;
	this->use_potek_kff = 0;
	this->kncrit=0.0;
}

// ------------------------ mean molecular weights ------------------------

double Eos::Yi(void)
// inverse mean molecular weight per ion
{
	if (set_Yi > 0.0) return set_Yi;
	double sum=0.0;
	for (int i=1; i<=this->ns; i++) sum+=this->X[i]/this->A[i];
	return sum;
}

double Eos::Ye(void)
// inverse mean molecular weight per electron
{
	if (set_Ye > 0.0) return set_Ye;
	double sum=0.0;
	for (int i=1; i<=this->ns; i++) sum+=this->X[i]*this->Z[i]/this->A[i];
	return sum;
}

double Eos::YZ2(void)
// sum of (X/A)*Z^2
{
	if (this->set_YZ2 > 0.0) return this->set_YZ2;
	double sum=0.0;
	for (int i=1; i<=this->ns; i++) sum+=this->X[i]*this->Z[i]*this->Z[i]/this->A[i];
	return sum;
}

// ------------------------ equation of state ------------------------------

double Eos::pe(void)
  // Calculates the (cgs) electron gas pressure as given by the
  // semi-analytic formula of Paczynski (1983) ApJ 267 315.
{
	double rY, pednr, pedr, pend, ped;
	rY=this->rho*Ye();
	pednr=9.91e-2*pow(rY,5.0/3.0);
	pedr=1.231e1*pow(rY,4.0/3.0);
	ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
	pend=8.254e-7*1e8*this->T8*rY;
	return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}

double Eos::ptot(void)
  // Calculates the total pressure, that is the sum of electron +
  // ion + radiation pressure. Coulomb corrections included if
  // species 1 has Z>2 i.e. not hydrogen or helium, this is a trick to only
  // apply Coulomb corrections in the ocean
{
	double P;

	if (this->use_potek_eos) {
		double dummy;
		potek_eos(&P,&dummy,&dummy);
	} else {
		
  double f;
  if (this->Z[1]>2.0) f=1.0+this->Uex()/3.0; else f=1.0;

  P=8.254e15*this->rho*this->T8*Yi()*f;   // ions
  P+=pe();                                // electrons  
  P+=RADa*1e32*pow(this->T8,4)/3.0;         // radiation
  //P+=2.521967e17*pow(this->T8,4);         // radiation
  if (this->Yn > 0.0) {                   // neutrons
    double k, EFn;
    // Assume the neutrons are non-relativistic
    // EFn=1.42*pow(1e-12*rho*Yn,2.0/3.0);
    // but we use a fit for EFn which includes interactions
    // comes from Mackie and Baym
    k=0.207*pow(1e-12*this->rho*this->Yn, 1.0/3.0);
    EFn=1.730*k+25.05*k*k-30.47*k*k*k+17.42*k*k*k*k;  // in MeV
    P+=0.4*EFn*1.6e-6*this->rho*this->Yn/1.66e-24;
  }
	}
  	return P;
}


double Eos::Utot(void)
  // internal energy density
{
	double f, r, tau;
	tau=this->T8/59.4;
	r=this->FermiI(2,this->T8,this->Chabrier_EF())/this->FermiI(1,this->T8,this->Chabrier_EF());
	if (this->Z[1]>1.0) f=1.5+this->Uex(); else f=1.5;
	return 8.254e15*this->rho*this->T8*Yi()*f + 1.5*this->pe()*(1.0+tau*r)/(1.0+0.5*tau*r) +
		3.0*2.521967e17*pow(this->T8,4);
}

double Eos::FermiI(int k, double T8, double EF)
  // fitting formula for the generalized Fermi integrals
  // from Chabrier and Potekhin 1998
  // nu=k+1/2, with k=0,1, or 2
{
  double c[3][5]={{0.37045057, 0.41258437, 9.777982e-2, 5.3734153e-3, 3.8746281e-5},
		  {0.39603109, 0.69468795, 0.22322760, 1.5262934e-2, 1.3081939e-4},
		  {0.76934619, 1.7891437, 0.70754974, 5.6755672e-2, 5.5571480e-4}};
  double e[3][5]={{0.43139881, 1.7597537, 4.1044654, 7.7467038, 13.457678},
		  {0.81763176, 2.4723339, 5.1160061, 9.0441465, 15.049882},
		  {1.2558461, 3.2070406, 6.1239082, 10.316126, 16.597079}};
  double x[5]={7.265351e-2, 0.2694608, 0.533122, 0.7868801, 0.9569313};
  double z[5]={0.26356032, 1.4134031, 3.5964258, 7.0858100, 12.640801};
  double h[5]={3.818735e-2, 0.1256732, 0.1986308, 0.1976334, 0.1065420};
  double v[5]={0.29505869, 0.32064856, 7.3915570e-2, 3.6087389e-3, 2.3369894e-5};
  int i;
  double I=0.0, F, R, chi, tau;

  tau=T8/59.4;
  chi=EF/(8.625*T8);

  R=sqrt(chi*(1.0+0.5*chi*tau));

  if (chi*tau > 0.01) {
    F=(chi+1.0/tau)*0.5*R-log(1+tau*chi+sqrt(2.0*tau)*R)/pow(2*tau,1.5);
    if (k>0) F=(2*pow(R,3.0)/3.0-F)/tau;
    if (k>1) F=(2*chi*pow(R,3.0)-5.0*F)/(4.0*tau);
  } else {
    F=pow(chi,1.0*k+1.5)/(1.0*k+1.5);
  }

  if (chi <= 0.6) {
    for (i=0; i<5; i++) {
      I+=c[k][i]*sqrt(1+0.5*e[k][i]*tau)/(exp(-e[k][i])+exp(-chi));
    }
  }
  if (chi > 0.6 && chi < 14.0) {
    for (i=0; i<5; i++) {
      I+=h[i]*pow(x[i],1.0*k)*pow(chi, 1.0*k+1.5)*sqrt(1+0.5*chi*x[i]*tau)/
	  (1+exp(chi*x[i]-chi));
      I+=v[i]*pow(z[i]+chi, 1.0*k+0.5)*sqrt(1+0.5*(z[i]+chi)*tau);
    }
  }
  if (chi >= 14.0) {
    I=F+(PI*PI/6.0)*pow(chi, 1.0*k)*(1.0*k+0.5+0.5*(k+1)*chi*tau)/R;
  }

  return I;
}


double Eos::Fermi_Inv_1_2(double F)
{
  double AN=0.5, RN, DEN, INV, FF;
  int i;
  int M1=2, K1=2, M2=2, K2=2;
  double A1[4]={0.0, 4.4593646e1, 1.1288764e1, 1.0};
  double B1[4]={0.0, 3.9519346e1, -5.7517464, 2.6594291e-1};
  double A2[4]={0.0, 3.4873722e1, -2.6922515e1, 1.0};
  double B2[4]={0.0, 2.6612832e1, -2.0452930e1, 1.1808945e1};
  
  if (F < 4.0) {
    RN=F+A1[M1];
    for (i=M1-1; i>=1; i--) RN=RN*F+A1[i];
    DEN=B1[K1+1];
    for (i=K1; i>=1; i--) DEN=DEN*F+B1[i];
    INV=log(F*RN/DEN);
  } else {
    FF=1.0/pow(F,1.0/(1.0+AN));
    RN=FF+A2[M2];
    for (i=M2-1; i>=1; i--) RN=RN*FF+A2[i];
    DEN=B2[K2+1];
    for (i=K2; i>=1; i--) DEN=DEN*FF+B2[i];
    INV=RN/(DEN*FF);
  }
  
  return INV;
}





double Eos::Chabrier_EF(void)
  // Calculates the Fermi energy in keV including the
  // rest mass using the fit of Chabrier and Potekhin,
  // (1998) Phys Rev E, 58, 4941
  // It uses Antia (1993) to evaluate the inverse Fermi integral
  //
  // rY is (rho/mu_e) in g/cm^3 ;  T is the temperature in K
{
  double EFnr, kT, x, tau, theta;
  double q1,q2,q3, et,etu, corr,F, mc2, rY, T;

  T=this->T8*1e8;
  rY=this->rho*Ye();

  // Electron rest mass in keV
  mc2=510.999;

  // Find kT, x=p_F/m_e c, tau=kT/me and theta=T/T_F
  kT=8.617347*T*1e-8;
  x=1.007e-2*pow(rY,1.0/3.0);
  tau=kT/mc2;
  theta=tau/(sqrt(1.0+x*x)-1.0);

  // Calculate the non-relativistic guess
  F=2.0*pow(theta,-1.5)/3.0;
  EFnr=kT*Fermi_Inv_1_2(F);

  // These functions are defined in CP eq. 24
  if (theta > 69.0) {
    et=1e30; etu=1e-30;
  } else {
    et=exp(theta); etu=1.0/et;
  }
  q1=1.5/(et-1.0);
  q2=12.0+8.0/pow(theta,1.5);
  q3=1.366-(etu+1.612*et)/(6.192*pow(theta,0.0944)*etu+
			   5.535*pow(theta,0.698)*et);

  // This is the correction to the non-relativistic EF
  corr=(1+q1*sqrt(tau)+q2*q3*tau)/(1+q2*tau);
  corr*=tau/(1+(tau/(2*theta)));
  corr=1.5*log(1+corr);

  // return E_F including the rest mass
  return mc2+EFnr-kT*corr;
}





void Eos::set_composition_by_density(void)
  // works out the composition at density rho according to the class variable 'accr'
  // accr=1 or 2 accreted crust; accr=0 equilibrium crust
{
	// accreted matter composition from Haensel & Zdunik (1990)
	double Acell[19]={56.0,56.0,56.0,56.0,56.0,56.0,56.0,56.0,112.0,112.0,112.0,112.0,112.0,224.0,224.0,224.0,224.0,448.0,448.0};
	double Aa[19]={56.0,56.0,56.0,56.0,56.0,52.0,46.0,40.0,68.0,62.0,56.0,50.0,44.0,66.0,60.0,54.0,48.0,96.0,88.0};
	double Za[19]={26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,12.0,18.0,16.0,14.0,12.0,24.0,22.0};
	double rhomaxa[19]={1.494e9,1.1145e10,7.848e10,2.496e11,6.110e11,9.075e11,1.131e12,1.455e12,1.766e12,2.134e12,2.634e12,3.338e12,4.379e12,5.665e12,7.041e12,8.980e12,1.127e13,1.137e13,1.253e13};
	// accreted matter composition from Haensel & Zdunik (2003)
	double Acell2[29]={106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,212.0,212.0,212.0,212.0,424.0,424.0,424.0,424.0,424.0,848.0,848.0,848.0,848.0};
	double Aa2[29]={106.0,106.0,106.0,106.0,106.0,106.0,106.0,92.0,86.0,80.0,74.0,68.0,62.0,56.0,50.0,42.0,72.0,66.0,60.0,54.0,92.0,86.0,80.0,74.0,68.0,124.0,120.0,118.0,116.0};
	double Za2[29]={44.0,42.0,40.0,38.0,36.0,34.0,32.0,28.0,26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,24.0,22.0,20.0,18.0,16.0,28.0,26.0,24.0,22.0};
  	double rhomaxa2[29]={3.517e8,5.621e9,2.413e10,6.639e10,1.455e11,2.774e11,4.811e11,7.785e11,8.989e11,1.032e12,1.197e12,1.403e12,1.668e12,2.016e12,2.488e12,3.153e12,3.472e12,4.399e12,5.355e12,6.655e12,8.487e12,9.242e12,1.096e13,1.317e13,1.609e13,2.003e13,2.520e13,3.044e14,3.844e13};
	// Haensel & Pichon (1994) for cold catalysed matter
	// from Table 1, missing off the last element
	double Ab[13]={56.0,62.0,64.0,66.0,86.0,84.0,82.0,80.0,78.0,126.0,124.0,122.0,120.0};
	double Zb[13]={26.0,28.0,28.0,28.0,36.0,34.0,32.0,30.0,28.0,44.0,42.0,40.0,38.0};
	double rhomaxb[13]={7.96e6,2.71e8,1.30e9,1.48e9,3.12e9,1.10e10,2.80e10,5.44e10,9.64e10,1.29e11,1.88e11,2.67e11,3.79e11};
	// continued using Douchin & Haensel (2001) Table 1 and 2  (nc is nb in units of fm-3; xc is the neutron fraction)
  	double Ac[43]={130.076,135.750,139.956,141.564,142.161,142.562,143.530,144.490,145.444,146.398,147.351,148.306,149.263,151.184,154.094,156.055,159.030,162.051,166.150,170.333,175.678,181.144,187.838,195.775,202.614,211.641,220.400,224.660,229.922,235.253,240.924,245.999,253.566,261.185,270.963,283.993,302.074,328.489,357.685,401.652,476.253,566.654,615.840};
	double Zc[43]={42.198,42.698,43.019,43.106,43.140,43.163,43.215,43.265,43.313,43.359,43.404,43.447,43.490,43.571,43.685,43.755,43.851,43.935,44.030,44.101,44.155,44.164,44.108,43.939,43.691,43.198,42.506,42.089,41.507,40.876,40.219,39.699,39.094,38.686,38.393,38.281,38.458,39.116,40.154,42.051,45.719,50.492,53.162};
	double nc[43]={1.2126e-4,1.6241e-4,1.9772e-4,2.0905e-4,2.2059e-4,2.3114e-4,2.6426e-4,3.0533e-4,3.5331e-4,4.0764e-4,4.6800e-4,5.3414e-4,6.0594e-4,7.6608e-4,1.0471e-3,1.2616e-3,1.6246e-3,2.0384e-3,2.6726e-3,3.4064e-3,4.4746e-3,5.7260e-3,7.4963e-3,9.9795e-3,1.2513e-2,1.6547e-2,2.1405e-2,2.4157e-2,2.7894e-2,3.1941e-2,3.6264e-2,3.9888e-2,4.4578e-2,4.8425e-2,5.2327e-2,5.6264e-2,6.0219e-2,6.4183e-2,6.7163e-2,7.0154e-2,7.3174e-2,7.5226e-2,7.5959e-2};
	double xc[43]={0.0000,0.0000,0.0000,0.0000,0.0247,0.0513,0.1299,0.2107,0.2853,0.3512,0.4082,0.4573,0.4994,0.5669,0.6384,0.6727,0.7111,0.7389,0.7652,0.7836,0.7994,0.8099,0.8179,0.8231,0.8250,0.8249,0.8222,0.8200,0.8164,0.8116,0.8055,0.7994,0.7900,0.7806,0.7693,0.7553,0.7381,0.7163,0.6958,0.6699,0.6354,0.6038,0.5898};

	double Z,A;
  	int i;
	switch (this->accr) {
		case 2:   // accreted composition (A=106)
    		i=0; while (this->rho > rhomaxa2[i] && i<29) i++;
    		if (i==29) i=28; // higher density than HZ's table, set it to the last value
			A=Aa2[i]; Z=Za2[i];
       		this->Yn=(Acell2[i]-A)/Acell2[i];
			break;
		case 1: // accreted composition (A=56)
    		i=0; while (this->rho > rhomaxa[i] && i<18) i++;
			A=Aa[i]; Z=Za[i];
			this->Yn=(Acell[i]-A)/Acell[i];
			break;
		default: // cold matter composition
      		if (this->rho < rhomaxb[11]) {
				i=0; while (this->rho > rhomaxb[i] && i<12) i++;
				A=Ab[i]; Z=Zb[i]; this->Yn=0.0;
      		} else {
				i=0; while (this->rho > nc[i]*1.66e15 && i<42) i++;
				A=Ac[i-1]+(Ac[i]-Ac[i-1])*(this->rho-1.66e15*nc[i-1])/(1.66e15*(nc[i]-nc[i-1]));
				Z=Zc[i-1]+(Zc[i]-Zc[i-1])*(this->rho-1.66e15*nc[i-1])/(1.66e15*(nc[i]-nc[i-1]));
				this->Yn=xc[i-1]+(xc[i]-xc[i-1])*(this->rho-1.66e15*nc[i-1])/(1.66e15*(nc[i]-nc[i-1]));
      		}
	}

	// We set A[1] to be Acell. This means that Yi = 1/Acell, correctly giving the ion density
  	this->X[1]=1.0; this->Z[1]=Z; this->A[1]=A/(1.0-this->Yn);
	this->set_Ye=(1.0-this->Yn)*Z/A;
}



void Eos::set_composition_by_pressure(void)
  // works out the composition at density rho according to the class variable 'accr'
  // accr=1 or 2 accreted crust; accr=0 equilibrium crust
{
	// accreted matter composition from Haensel & Zdunik (1990)
	double Acell[19]={56.0,56.0,56.0,56.0,56.0,56.0,56.0,56.0,112.0,112.0,112.0,112.0,112.0,224.0,224.0,224.0,224.0,448.0,448.0};
	double Aa[19]={56.0,56.0,56.0,56.0,56.0,52.0,46.0,40.0,68.0,62.0,56.0,50.0,44.0,66.0,60.0,54.0,48.0,96.0,88.0};
	double Za[19]={26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,12.0,18.0,16.0,14.0,12.0,24.0,22.0};
	double Pmaxa[19] = {7.235e26,9.569e27,1.152e29,4.747e29,1.361e30,1.980e30,2.253e30,2.637e30,2.771e30,3.216e30,3.825e30,4.699e30,6.043e30,7.233e30,9.238e30,1.228e31,1.602e31,1.613e31,1e33};
	// accreted matter composition from Haensel & Zdunik (2003)
	double Acell2[30]={106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,106.0,212.0,212.0,212.0,212.0,424.0,424.0,424.0,424.0,424.0,848.0,848.0,848.0,848.0,848.0};
	double Aa2[30]={106.0,106.0,106.0,106.0,106.0,106.0,106.0,92.0,86.0,80.0,74.0,68.0,62.0,56.0,50.0,42.0,72.0,66.0,60.0,54.0,92.0,86.0,80.0,74.0,68.0,124.0,120.0,118.0,116.0,116.0};
	double Za2[30]={44.0,42.0,40.0,38.0,36.0,34.0,32.0,28.0,26.0,24.0,22.0,20.0,18.0,16.0,14.0,12.0,20.0,18.0,16.0,14.0,24.0,22.0,20.0,18.0,16.0,28.0,26.0,24.0,22.0,22.0};
	double Pmaxa2[30]={9.235e25,3.603e27,2.372e28,8.581e28,2.283e29,5.025e29,9.713e29,1.703e30,1.748e30,1.924e30,2.135e30,2.394e30,2.720e30,3.145e30,3.723e30,4.549e30,4.624e30,5.584e30,6.883e30,8.749e30,1.157e31,1.234e31,1.528e31,1.933e31,2.510e31,3.363e31,4.588e31,5.994e31,8.408e31,1e34};
	// Haensel & Pichon (1994) for cold catalysed matter
	// from Table 1, missing off the last element
	double Ab[13]={56.0,62.0,64.0,66.0,86.0,84.0,82.0,80.0,78.0,126.0,124.0,122.0,120.0};
	double Zb[13]={26.0,28.0,28.0,28.0,36.0,34.0,32.0,30.0,28.0,44.0,42.0,40.0,38.0};
	double Pmaxb[13]={5.26021e+23, 6.94463e+25, 5.4887e+26, 6.26676e+26, 1.65627e+27, 8.54559e+27, 2.83828e+28, 6.542e+28, 1.32675e+29, 1.85549e+29, 2.9496e+29, 4.27004e+29, 6.09914e+29};
	// continued using Douchin & Haensel (2001) Table 1 and 2  (nc is nb in units of fm-3; xc is the neutron fraction)
  	double Ac[43]={130.076,135.750,139.956,141.564,142.161,142.562,143.530,144.490,145.444,146.398,147.351,148.306,149.263,151.184,154.094,156.055,159.030,162.051,166.150,170.333,175.678,181.144,187.838,195.775,202.614,211.641,220.400,224.660,229.922,235.253,240.924,245.999,253.566,261.185,270.963,283.993,302.074,328.489,357.685,401.652,476.253,566.654,615.840};
	double Zc[43]={42.198,42.698,43.019,43.106,43.140,43.163,43.215,43.265,43.313,43.359,43.404,43.447,43.490,43.571,43.685,43.755,43.851,43.935,44.030,44.101,44.155,44.164,44.108,43.939,43.691,43.198,42.506,42.089,41.507,40.876,40.219,39.699,39.094,38.686,38.393,38.281,38.458,39.116,40.154,42.051,45.719,50.492,53.162};
	double xc[43]={0.0000,0.0000,0.0000,0.0000,0.0247,0.0513,0.1299,0.2107,0.2853,0.3512,0.4082,0.4573,0.4994,0.5669,0.6384,0.6727,0.7111,0.7389,0.7652,0.7836,0.7994,0.8099,0.8179,0.8231,0.8250,0.8249,0.8222,0.8200,0.8164,0.8116,0.8055,0.7994,0.7900,0.7806,0.7693,0.7553,0.7381,0.7163,0.6958,0.6699,0.6354,0.6038,0.5898};
	double Pmaxc[43]={3.09981e+29, 4.31891e+29, 5.44304e+29, 5.78936e+29, 5.99198e+29, 6.1361e+29, 6.53641e+29, 7.00062e+29, 7.53879e+29, 8.15861e+29, 8.87021e+29, 9.67731e+29, 1.05913e+30, 1.27633e+30, 1.69967e+30, 2.05599e+30, 2.71729e+30, 3.5506e+30, 4.97197e+30, 6.8044e+30, 9.78126e+30, 1.36685e+31, 1.97939e+31, 2.94326e+31, 4.03465e+31, 5.96621e+31, 8.56788e+31, 1.01615e+32, 1.24515e+32, 1.50802e+32, 1.80384e+32, 2.06167e+32, 2.40494e+32, 2.69024e+32, 2.97801e+32, 3.25884e+32, 3.52151e+32, 3.74878e+32, 3.88157e+32, 3.96266e+32, 3.96055e+32, 3.87876e+32, 3.82389e+32};

	double Z,A;
  	int i;
	switch (this->accr) {
		case 2: {  // accreted composition (A=106)
    		i=0; while (this->P > Pmaxa2[i] && i<30) i++;
			A=Aa2[i]; Z=Za2[i];
       		this->Yn=(Acell2[i]-A)/Acell2[i];
			} break;
		case 1: { // accreted composition (A=56)
    		i=0; while (this->P > Pmaxa[i] && i<18) i++;
			A=Aa[i]; Z=Za[i];
			this->Yn=(Acell[i]-A)/Acell[i];
			} break;
		default: // cold matter composition
      		if (this->P < Pmaxb[11]) {
				i=0; while (this->P > Pmaxb[i] && i<12) i++;
				A=Ab[i]; Z=Zb[i]; this->Yn=0.0;
      		} else {
				i=0; while (this->P > Pmaxc[i] && i<42) i++;
				if (i==42) {
					A=Ac[i];
					Z=Zc[i];
					this->Yn=xc[i];
				} else {
					A=Ac[i-1]+(Ac[i]-Ac[i-1])*(this->P-Pmaxc[i-1])/(Pmaxc[i]-Pmaxc[i-1]);
					Z=Zc[i-1]+(Zc[i]-Zc[i-1])*(this->P-Pmaxc[i-1])/(Pmaxc[i]-Pmaxc[i-1]);
					this->Yn=xc[i-1]+(xc[i]-xc[i-1])*(this->P-Pmaxc[i-1])/(Pmaxc[i]-Pmaxc[i-1]);
				}
      		}
	}

	// We set A[1] to be Acell. This means that Yi = 1/Acell, correctly giving the ion density
  	this->X[1]=1.0; this->Z[1]=Z; this->A[1]=A/(1.0-this->Yn);
	this->set_Ye=(1.0-this->Yn)*Z/A;
}





// ----------------------- thermodynamics --------------------------------

double Eos::chi(double *x)
{
	double x1, x2, p1, p2;
	x1=*x; p1=ptot();
	x2=*x=1.001*x1; p2=ptot();
	*x=x1;
	return (log(p2)-log(p1))/(log(x2)-log(x1));
}

double Eos::CP(void)
// Calculates specific heat at constant pressure
{
	double chirho=chi(&this->rho);
	double chiT=chi(&this->T8);
	return CV()+chiT*chiT*ptot()/(this->rho*this->T8*1e8*chirho);
}

double Eos::Gamma1(void)
{
	double chirho=chi(&this->rho);
	double chiT=chi(&this->T8);
	return chirho+chiT*chiT*ptot()/(CV()*this->rho*1e8*this->T8);
}  

double Eos::del_ad(void)
  // calculates dlnT/dlnp at constant entropy
{
	double chirho=chi(&this->rho);
	double chiT=chi(&this->T8);
	double gam1=chirho+chiT*chiT*ptot()/(CV()*this->rho*1e8*this->T8);
	return ptot()*chiT/(CV()*gam1*this->rho*1e8*this->T8);
}

double Eos::CV(void)
// Calculates the specific heat at constant volume (density)
{
	if (this->use_potek_eos && !(this->Yn>0.0)) {   // use Potekhin's eos only before neutron drip
		double dummy, dummy2, dummy3;
		potek_eos(&dummy,&dummy2,&dummy3);
		this->cvion = dummy3;
		this->cve=dummy2;
	} else {
		
		// IONS		
		double gg=this->gamma();

		if (gg < this->gamma_melt) {  // liquid

			// alpha is from Potekhin & Chabrier 2000
			double a1,a2,a3,b1,b2,b3,b4;
			a1=-0.9070; a2=0.62954; a3=-0.5*sqrt(3.0)-a1/sqrt(a2);
			b1=4.56e-3; b2=211.6; b3=-1.0e-4; b4=4.62e-3;
			double alpha=0.5*pow(gg,1.5)*(a3*(gg-1.0)/pow(gg+1.0,2.0)-a1*a2/pow(gg+a2,1.5))
				+pow(gg,2.0)*(b3*(pow(gg,2.0)-b4)/pow(pow(gg,2.0)+b4,2.0)-b1*b2/pow(gg+b2,2.0));
			this->cvion=8.3144e7*(1.5+alpha)*this->Yi();
			cv_alpha=alpha;

  		} else {  // solid    -- This implements equation (5) of Chabrier 1993 
	
			// In the next line, I am putting the mass of the nucleus  A = Acell(1-Yn) = A[1](1-Yn), ie. no entrainment
			double eta=7.76e-5*this->Z[1]*sqrt(this->Yi()*this->rho/(this->A[1]*(1.0-this->Yn)))/this->T8;
			double alphaeta=0.399*eta;
			double gameta=0.899*eta;
			double dd,dd1,dd2;
			double x=alphaeta;
			dd1=pow(3.141592654,4.0)/(5.0*pow(x,3.0));
			dd1-=3.0*exp(-x)*(6.0+x*(6.0+x*(3.0+x)))/pow(x,3.0);
			dd2=1.0-0.375*x+0.05*x*x;
			if (dd1 > dd2) dd=dd2; else dd=dd1;
			this->cvion=8.3144e7*this->Yi()*(8.0*dd-6*alphaeta/(exp(alphaeta)-1.0)+(pow(gameta,2.0)*exp(gameta)/pow(exp(gameta)-1.0,2.0)));	
			if (isnan(this->cvion)) this->cvion=0.0;

		}
  	
  		// ELECTRONS
		{ // modified version of Paczynksi's fit for cve
    		double dT,temp,p1,p2;
			temp=1.001*this->T8; dT=temp-this->T8;
       		p1=this->pemod(); this->T8+=dT; p2=this->pemod(); this->T8-=dT;
			this->cve=(1/((this->f()-1)*this->rho))*1e-8*(p2-p1)/dT;
		}
	}

  	// RADIATION
    this->cvrad=4.0*RADa*1e24*pow(this->T8,3)/this->rho;

  	// NEUTRONS
	this->cvneut=0.0;
  	if (this->Yn > 0.0) {
		double EFn;
		EFn=1.42*pow(1e-12*rho*Yn,2.0/3.0);
		cvneut=0.5*3.1415*3.1415*8.3144e7*this->Yn * 1.38e-16*1e8*this->T8/(EFn*1.6e-6);
		if (this->gap > 0) {
			double R00;
			double t, u;
			t = 1e8*this->T8/this->TC();
			if (t>1.0) t=1.0;
			u = sqrt(1.0-t)*(1.456 - 0.157/sqrt(t)+1.764/t);			
			R00 = pow(0.4186+sqrt(1.007*1.007 + pow(0.5010*u,2.0)),2.5) * exp(1.456 - sqrt(1.456*1.456+u*u));
			cvneut *= R00;
		}
	}

	return this->cvion+this->cve+this->cvrad+cvneut;
}


double Eos::pemod(void)
// Modified version of the electron pressure formula to use for heat capacity (Paczynski 1983 ApJ 267 315)
{
	double rY, pednr, pedr, pend, ped;
	rY=this->rho*this->Ye();
	// divide pednr and pedr by appropriate factors
	pednr=9.91e-2*pow(rY,5.0/3.0)/1.32;
	pedr=1.231e1*pow(rY,4.0/3.0)/0.822;
	ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
	pend=8.254e-7*1e8*this->T8*rY;
	return(1e14*sqrt(pow(ped,2)+pow(pend,2)));
}

double Eos::f(void)
// Calculates f = dln ped/dln rho, using the fitting formula given by Paczynski (1983)
{
	double rY, pednr, pedr, ped;
	rY=this->rho*Ye();
	pednr=9.91e-2*pow(rY,5.0/3.0);
	pedr=1.231e1*pow(rY,4.0/3.0);
	ped=1/sqrt((1/pow(pedr,2))+(1/pow(pednr,2)));
	return (5.0*pow(ped/pednr,2) + 4.0*pow(ped/pedr,2))/3.0;
}

// --------------------------- opacity ------------------------------------


double Eos::eps_nu(void)
// Calculates neutrino emissivity (erg/g/s)
// by plasma process from Schinder et al. 1987
// and includes neutrino synchtron from Bezchastnov et al. 1997 (relevant for high B)
{
	double a0, a1, a2, b1, b2, b3, c;
	double xi, xi2, xi3, la, la2, la3, g, K;

	Q6=0.0;   // Q6 is not included in the total emissivity, 
	// it can be used to return a comparison value

	// variables
	la=this->T8/59.302; la2=la*la; la3=la2*la;
	xi=pow(this->rho*this->Ye()*1e-9, 1.0/3.0)/la;
	xi2=xi*xi; xi3=xi2*xi;

	double xx = this->x();
//	xx = pow(this->rho*this->Ye()/9.7e5,1.0/3.0);

	// ------ 1. plasma ------
	// these coefficients valid for 10^8<T<10^11 K
	a0=2.146e-7; a1=7.814e-8; a2=1.653e-8;
	b1=2.581e-2; b2=1.734e-2; b3=6.990e-4;
	c=0.56457;

	K=pow(this->rho*this->Ye(),3.0);

	// formula from Schinder et al.
	Q1=K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));

	// Also check the formula from Yakovlev et al. (2001) eqs (23,38,40)
	// At high density, when the emissivity is suppressed the agreement is excellent
	// At low density, Yakovlev et al. are tens of percent to factor of 2 larger
	if (0) {
		double Qc=1.023e23;
		double xr = 100.9*pow(1e-12*this->rho*this->Ye(),1.0/3.0);
		double tr = 1.38e-8*this->T8/(9.11e-28*9e20);
		double fp = sqrt(4.0*(1.0/137.0)*pow(xr,3.0)/(3.0*PI*sqrt(1.0+xr*xr)))/tr;
		double Ip = pow(tr,9.0)*(16.23*pow(fp,6.0)+4.604*pow(fp,15.0/2.0))*exp(-fp);
		Q6 = Ip*Qc*0.9248*137.0/(96.0*pow(PI,4.0));
	}
	
	
	// ------ 2. pair ------
	// coefficients valid for 10^8 < T < 10^11 K
	a0=5.026e19; a1=1.745e20; a2=1.568e21;
	if (this->T8 < 100.0) {  // 10^8<T<10^10 K
		b1=9.383e-1; b2=-4.141e-1; b3=5.829e-2;
		c=5.5924;
	} else { // 10^10 < T < 10^11 K
		b1=1.2383; b2=-8.141e-1; b3=0.0;
		c=4.9924;
	}

	g=1.0-13.04*la2+133.5*la2*la2+1534*la2*la2*la2+918.6*la2*la2*la2*la2;
	K=g*exp(-2.0/la);

	double qpair;
	qpair=pow(10.7480*la2+0.3967*sqrt(la)+1.0050,-1.0)
		* pow(1.0 + this->rho*this->Ye()/(7.692e7*la3+9.715e6*sqrt(la)),-0.3);

	// formula from Schinder et al.
	Q2=(1.0+0.10437*qpair)*K*exp(-c*xi)*(a0+a1*xi+a2*xi2)/(xi3+(b1/la)+(b2/la2)+(b3/la3));

	
	// ------ 3. Brems ------
	
	//Q3=0.3229*this->rho*this->YZ2()*pow(this->T8,6.0);

	if (this->gamma() < this->gamma_melt) { // liquid
		// Here we implement equations (8) and (25) from Haensel et al. 96
		double L;
		double A, B, eta, t, Z;
		Z=this->Ye()/this->Yi();
		t=this->T8/(118.6*(sqrt(1.0+pow(xx,2.0))-1.0));

		// finite core radius : value depends on whether above or
		// below neutron drip
		if (this->Yn == 0.0) eta=0.16*pow(this->rho*1e-12,1.0/3.0);
		else eta=0.25*pow(this->rho*1e-12*this->Ye(),1.0/3.0);
		// Eq.(25) of Haensel et al. 1996 for the Coulomb logarithm L
		A=0.269+20.0*t+0.0168*Z+0.00121*eta-0.0356*Z*eta+0.0137*Z*Z*t+1.54*Z*t*eta;
		B=1.0+180.0*t*t+0.483*t*Z+20.0*t*Z*eta*eta+4.31e-5*Z*Z;
		L=A/pow(B,0.75);

		Q3=L*0.3229*this->rho*this->YZ2()*pow(this->T8,6.0);

	} else {  // solid

		// Fit from Kaminker et al. (1999), eq.(40)
		double r=log10(this->rho*1e-12);
		double t=log10(this->T8);
		double r0=2.8e14;

		Q3=11.204 + 7.304*t + 0.2976*r - 0.370*t*t + 0.188*t*r - 0.103*r*r
			+ 0.0547*t*t*r - 6.77*log10(1+0.228*this->rho/r0);

		if (this->accr) {
		// adjustment to approximate accreted matter
		// as described in Cumming et al. (2006)
			if (r < 0.1) Q3-=0.2;
			else {
				if (r < 1.0) Q3-=0.3;
				else Q3-=0.4;
			}
		}

		Q3=pow(10.0,Q3);
	}


	// ------ 4. Cooper pair emission in the crust ------
	Q4=0.0;
	if (0) {   // switch it off
		if (this->Yn > 0.0 && this->T8*1e8<TC()) {
			double kf=0.261*pow(1e-12*this->rho*this->Yn,1.0/3.0);
			double pf=kf*197.0;
			double tau=this->T8*1e8/this->TC();
			double u=sqrt(1.0-tau)*(1.456-(0.157/sqrt(tau))+(1.764/tau));
			double fac=0.602*u*u+0.5942*pow(u,4.0)+0.288*pow(u,6.0);
			fac*=sqrt(0.5547+sqrt(pow(0.4453,2.0)+0.01130*u*u));
			fac*=exp(-sqrt(4*u*u+pow(2.245,2.0))+2.245);

			Q4=3.0*1.17e21*(pf/940.0)*pow(this->T8*0.1,7.0);
			Q4*=fac;

			Q4*=this->Yn/(this->Yn+(1.0-this->Yn)/this->A[1]);  
			// ^^^need to multiply by fraction of space
			// occupied by free neutrons. I think this is n_n/n_particles
			// this factor is not in Dany's code! but it's close to 1 anyway
			// ie. most of the particles are neutrons
		}
	}

	
	// ------ 5. Neutrino synchrotron ------
	//  Bezchastnov et al. 1997
	Q5 = 0.0;
	if (1) {  // use this to turn off neutrino synchrotron

		// First regime B given by their eq. 7
		// this is the simple density independent part
		Q5 = 9.04e14*pow(this->B/1e13,2.0)*pow(this->T8/10.0,5.0);

		// Next need to include suppression factors from transition to A and C
		double TB,z,xi;
		TB = 1.34e9*this->B*1e-13/sqrt(1.0+xx*xx);
		z = 1e-8*TB/this->T8;
		xi = 1.5*z*pow(xx,3.0);  // this is TP/T
		
		double SAB, SBC;
		{
			double D1,D2;
			D1 = 1.0+0.4228*z+0.1014*z*z+0.006240*z*z*z;
			D2 = 1.0+0.4535*pow(z,2.0/3.0)+0.03008*z-0.05043*z*z+0.004314*z*z*z;
			SBC = exp(-0.5*z)*D1/D2;
		}
		{
			double Fm,Fp, y1,y2;
			double a1=2.036e-4, b1=7.405e-8, c1=3.675e-4;
			double a2=3.356e-3, b2=1.536e-5, c2=1.436e-2, d2=1.024e-5, e2=7.647e-8;
			double DD1 = 44.01, DD2 = 36.97, alpha1=3172.0, alpha2=172.2;
			y1 = pow(pow(1.0+alpha1*pow(xi,2.0/3.0),2.0/3.0)-1.0,1.5);
			y2 = pow(pow(1.0+alpha2*pow(xi,2.0/3.0),2.0/3.0)-1.0,1.5);
			Fp = DD1 * pow(1.0+c1*y1,2.0)/pow(1.0+a1*y1+b1*y1*y1,4.0);
			Fm = DD2 * (1.0+c2*y2+d2*y2*y2+e2*y2*y2*y2)/pow(1.0+a2*y2+b2*y2*y2,5.0);	
			SAB = (27.0*pow(xi,4.0)/(PI*PI*512.0*1.037))*(Fp-0.175*Fm/1.675);
		}
		
		//SAB=1.0; SBC=1.0;  // turn off regimes A and C
		Q5 *= SBC*SAB;
	}

	// return the summed emissivity (divide by rho to get per gram)
	return (Q1+Q2+Q3+Q4+Q5)/this->rho;
}



double Eos::TC(void)
// calculates the critical temp for neutron superfluidity in the crust
{
	// neutron k_F in fm^-1
	double k=0.261*pow(1e-12*this->rho*this->Yn,1.0/3.0);
	double Tcrit;
 	
	switch (this->gap) {

		case 1: {// SFB03
  			double k0[18]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.175,1.25, 1.3, 1.35, 1.4, 1.45};
  			double d0[18]={0.0, 0.09, 0.210, 0.360, 0.500, 0.610, 0.720, 0.790, 0.780,0.700, 0.580, 0.450, 0.280, 0.190, 0.100, 0.030, 0.0};
  			double d2[18]={2.915e1, -4.297, 6.040, -1.863, -4.59, 2.221, -4.296,-9.037, -7.555, -2.741, -5.480, -1.344e1, 1.656e1, -6.667,1.010e1, 1.426e1, 2.887e1};
  			if (k < k0[0]) return 0.0;
  			if (k > k0[16]) return 0.0;
	  		int i1=0;
  			int i2=16;
  			int i;
  			while ((i2-i1)>1) {
    			i=(i2+i1)/2;
    			if (k0[i] > k) {
      				i2=i;
    			} else {
      				i1=i;
    			}
  			}
  			double delk=k0[i2]-k0[i1];
  			double a=(k0[i2]-k)/delk;
  			double b=(k-k0[i1])/delk;
  			double t=a*d0[i1]+b*d0[i2]+((pow(a,3.0)-a)*d2[i1]+(pow(b,3.0)-b)*d2[i2])*(delk*delk)/6;
  			Tcrit=(t/1.76)*1.1604e10;
			} break;
		
		case 2: { // AWP II
  			double k0[18]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,1.3, 1.4, 1.5, 1.6, 1.7};
  			double d0[18]={0.0, 3.3e8, 1.18e9, 2.44e9, 4.20e9, 6.13e9, 7.91e9,9.10e9, 9.56e9, 9.03e9, 7.71e9, 5.93e9,4.15e9, 2.50e9, 1.12e9, 3.61e8, 0.0};
  			double d2[18]={7.58e10, 4.64e10, 5.06e10, -2.87e9, 1.17e11, -7.45e10, -5.28e10,-6.83e10, -1.12e11, -7.75e10, -5.19e10, 9.16e9, 1.53e10, 7.71e9,1.16e11, -3.87e10, 1.58e11};
  			if (k < k0[0]) return 0.0;
  			if (k > k0[16]) return 0.0;
  			int i1=0;
  			int i2=16;
  			int i;
  			while ((i2-i1)>1) {
    			i=(i2+i1)/2;
    			if (k0[i] > k) {
      				i2=i;
    			} else {
      				i1=i;
    			}
  			}
  			double delk=k0[i2]-k0[i1];
  			double a=(k0[i2]-k)/delk;
  			double b=(k-k0[i1])/delk;
			double t=a*d0[i1]+b*d0[i2]+((pow(a,3.0)-a)*d2[i1]+(pow(b,3.0)-b)*d2[i2])*(delk*delk)/6;
  			Tcrit=t;
			} break;
		
		case 3: { // Gaussian Tc
   			double kbar=1.3;
   			double width=0.5;
			Tcrit=6e9*exp(-pow(k-kbar,2.0)/(2.0*width*width));
			} break;
		
		case 4: { // "all or nothing" cutoff for k<kc
			double k=0.261*pow(1e-12*this->rho*this->Yn,1.0/3.0);		
			if (k < this->kncrit) Tcrit=0.0; else Tcrit=1e10;
			} break;
		
		case 5: { //"B1" from Reddy and Page (2011) book chapter Fig 4
			double rh[36]={11.8,11.8211,11.8632,11.9265,12.0056,12.0532,12.0904,12.1223,12.1543,12.1914,12.2233,12.2657,12.3133,12.3609,12.3981,12.6160,12.7330,12.8235,12.8977,12.9721,13.0620,13.1781,13.2728,13.3251,13.3773,13.4136,13.4499,13.4808,13.5223,13.5480,13.5998,13.6617,13.7706,13.8536,13.9266,14.0};
			double T[36]={0.0,0.0319625,0.0961741,0.192491,0.451774,0.809955,1.23378,1.72311,2.27779,2.73430,3.12559,3.51659,3.90745,4.26563,4.72213,7.69005,9.45157,10.8543,11.6693,12.6150,13.2661,13.7205,13.5872,13.0956,12.5059,11.8187,10.9680,10.0848,9.16864,8.25291,7.07503,5.34132,3.05077,1.28380,0.334101,0.0};
			double rr=log10(this->rho);
			if (rr <= rh[0] || rr >= rh[35]) return 0.0;
			int ii=0; 
			while (rr>rh[ii]) ii++;
			double Tc=T[ii-1] + T[ii]*(rr-rh[ii-1])/(rh[ii]-rh[ii-1]);
			Tcrit=1e9*Tc;
			} break;
		
		case 6: { //BCS
			double Tc;
			if (k >0.0) {
				double xx = PI / (2.0*18.5*k);	
				Tc = (8.0/137.0) * exp(-xx) * pow(1.05e-27*1e13*k,2.0)/(2.0*1.67e-24);
				Tc/=1.38e-16;
			} else Tc=0.0;
			Tcrit=Tc/1.76;
			} break;
		
		default:   // Tc=0 everywhere:  no superfluidity
			Tcrit=0.0;
	}

	return Tcrit;
}




double Eos::opac(void)
  // Calculates the opacity
{
	// Electron scattering opacity from Paczynski
	this->kes=(0.4*Ye())/((1+2.7e11*this->rho*pow(1e8*this->T8,-2.0))*(1+pow(this->T8/4.5,0.86)));

	// Fermi energy
	double eta, ef=Chabrier_EF();
	if (ef == 0) {
	  eta=Fermi_Inv_1_2(1.105e-4*this->rho*Ye()/pow(this->T8,1.5));
	  ef=eta*8.617*this->T8;
	} else {
	  eta=(ef-me)/(8.617*this->T8);
	}

	// Free-Free opacity
	this->kff=7.53e-6*this->rho*Ye()/pow(this->T8, 3.5);
	kgaunt=0.0;
	for (int i=1; i<=this->ns; i++) kgaunt+=this->Z[i]*this->Z[i]*this->X[i]*gff(this->Z[i],eta)/this->A[i];
	this->kff*=kgaunt;
	if (this->use_potek_kff) {
		// Free-free opacity from Potekhin's magnetized envelopes
		double TRy = 100.0*this->T8/(0.15789*this->Z[1]*this->Z[1]);
		double c7 = 108.8 + 77.6*pow(TRy,0.834);
		c7 /= 1.0+0.502*pow(TRy,0.355)+0.245*pow(TRy,0.834);
		this->kff = this->kes * 2e4 * pow(this->Z[1],2.0) * this->rho /(c7 * this->A[1] * pow(100.0*this->T8,3.5));
	}
  	
	// total radiative opacity
	this->kappa_rad=this->kff+this->kes;

	// "non-additivity" factor from Potekhin et al. (2001) eqs 19-20
	double f=this->kff/this->kappa_rad;
	double TRy = 100.0*this->T8/(0.15789*this->YZ2()/this->Yi());
	double A = 1.0 + (1.097+0.777*TRy)*pow(f,0.617)*pow(1.0-f,0.77)/(1.0+0.536*TRy);
	this->kappa_rad*=A;

	if (this->B > 0.0) {
		// magnetic correction, Potekhin et al. 2001 eqs 21-23
		double xr;
		if (this->rho < 7.09e3*pow(1e-12*this->B,1.5)/this->Ye()) {
			xr = 2.96e-5*this->rho*this->Ye()*1e12/this->B;			
		} else {
			xr = this->x();
		}		
		double TB8 = 1.343*1e-12*this->B/sqrt(1.0+xr*xr);
		double u = TB8/(2.0*this->T8);
		double AA1,AA2,AA3;
		double a1=0.0949, a2=0.1619, a3=0.2587, b1=0.0610,
			b2=0.1400, b3=0.1941, c1=0.09, c2=0.0993, c3=0.0533;
		AA1=a1-b1*pow(f,c1);
		AA2=a2-b2*pow(f,c2);
		AA3=a3-b3*pow(f,c3);
		double corr = 1.0  + u*u*(AA1*u+pow(AA2*u,2.0))/(1.0+AA3*u*u);
		this->kappa_rad /= corr;
	}
	
	// Correction for plasma frequency from Potekhin et al. (2003) ApJ 594,404
	kappa_rad *= exp(0.005*log(1.0 + 1.5*sqrt(this->rho*1e-6*this->Ye())*(28.8/(8.625*this->T8))));	
	
  	// Conduction
	double KK;
	if (use_potek_cond) KK = potek_cond(); else KK = K_cond(ef);
  	this->kcond=3.024e20*pow(this->T8,3)/(KK*this->rho);
 
  	// Add up opacities in parallel
	return 1.0/((1.0/this->kcond)+(1.0/this->kappa_rad));
}


double Eos::gff(double Z1, double eta)
// Calculates the free-free Gaunt factor for element with
// charge Z1 using a fitting formula described in Schatz et al. (1999)
{
	double gaunt, x, rY, T8_32, gam;
	rY=this->rho*Ye();
	T8_32=pow(this->T8, 1.5);
	if (eta < 100.0) x=log(1.0+exp(eta)); else x=eta; // make sure it doesn't freak out for extremely large eta
	gaunt=1.16*8.02e3*x*T8_32/rY;  // normalisation and degeneracy piece
	x=pow(1+x,2.0/3.0);
	gam=sqrt(1.58e-3/this->T8)*Z1;
	gaunt*=(1.0-exp(-2*PI*gam/sqrt(x+10.0)))/(1.0-exp(-2*PI*gam/sqrt(x))); // Elwert factor
	gaunt*=1.0+pow(this->T8/7.7, 1.5);  // relativistic piece
	return gaunt;
}


double Eos::Uex(void)
  // Coulomb correction Uex/kT
{
	double u,g,g2,g3,g14;
	g=this->gamma(); g2=g*g; g3=g2*g; g14=pow(g,0.25);
	if (g < this->gamma_melt) {
		u=-0.89813*g+0.98686*g14-0.91095+0.25098/g14;
	} else {
		u=-0.89593*g+1.5+9.65/g+840/g2+1.101e5/g3;
	}
	return u;
}


double Eos::Fep(int flag)
  // "Coulomb log" for electron-phonon scattering 
  // (Baiko & Yakovlev 1995,1996)
  // if flag=0 electrical conductivity; flag>0 thermal conductivity
{
  double R0, R1, R2, G0, G2, t, u2, u1, s, F, K0;
  double alpha, alpha0, a0, a2, x, beta, AA, ZZ;
  double P0,g, K2,P2, c1, c2;

  // constants -- use values for bcc crystal
  a0=0.0174; a2=0.0118; u2=13.0; u1=2.8;

  AA=this->A[1]; ZZ=this->Z[1];
  x=this->x(); beta=x/sqrt(1+x*x);
  t=0.804*this->T8*(0.5*AA/ZZ)/sqrt(1e-9*this->rho);
  s=pow(4*ZZ,-2.0/3.0)+2.323e-3/beta;
  alpha0=1.683*sqrt(x/(AA*ZZ));
  alpha=alpha0*(0.5*u1*exp(-9.1*t)+t*u2);
  //alpha=1e-6;  //  small alpha is the Yakovlev & Urpin limit
  
  G0=u2*t/sqrt(t*t+a0); 
  R0=(exp(-alpha*s)-exp(-alpha))/alpha;
  R1=2*(exp(-alpha*s)*(1+alpha*s)-exp(-alpha)*(1+alpha))/
    (alpha*alpha);
  K0=2*R0-beta*beta*R1;
  
  // correction for finite nuclear size
  if (this->rho < 4e11) g=0.16*pow(this->rho*1e-12,1.0/3.0);
  else g=0.25*pow(this->rho*1e-12*this->Ye(),1.0/3.0);
  //  g=0.0;  switch off finite size effects
  P0=4.787-0.0346*ZZ;
  R2=(exp(-alpha*s)*(alpha*alpha*s*s+2*alpha*s+2)-exp(-alpha)*(alpha*alpha+2*alpha+2))/(alpha*alpha*alpha);
  c1=pow(1.0+pow(18.0*ZZ*PI,2.0/3.0)*g*g*(0.5*R1-beta*beta*R2)/(2.5*K0*P0),-P0);
  
  F=G0*K0*c1;
  
  if (flag > 0) { // thermal conductivity so add an extra piece
    P2=2.729-0.0204*ZZ;
    R2=this->Eep(alpha*s)-this->Eep(alpha);
    G2=t/(PI*PI*pow(t*t+a2,1.5));
    K2=0.5*R2-0.5*beta*beta*R0;
    // correction for finite nuclear size
    c2=pow(1.0+pow(18.0*PI*ZZ,2.0/3.0)*g*g*0.5*K0/(10.0*K2*P2),-P2);
    F+=G2*(3*K2-0.5*K0)*c2;
  }
  return F;
}

double Eos::Eep(double q)
  // used by Fep() to calculated thermal conductivity piece
  // Baiko & Yakovlev 1995
{
  double q2,q3,q4,qu;
  q2=q*q; q3=q2*q; q4=q3*q; qu=1.0/q;
  return exp(-q4/(q3+0.1397))*(log(1+qu)-0.5772/(1+2.2757*q2));
}


double Eos::lamei(int n)
// n=1 for thermal n=0 for electrical
{

  double L1, L2, s, w, I1, I2;


  if (this->T8 < 0.0) return 1.0;

if (this->rho > 3e14) return 1.0;

  if (0) {
  /* Yakovlev & Urpin */

  double x,x1,x2,lam;
  x1=this->x();
  x2=0.22*sqrt(this->T8);
  if (x1>x2) x=x1; else x=x2;
  lam=127*x*sqrt((3.0/this->gamma())+1.5)/
    pow(this->rho*this->Yi(),1.0/3.0);
  return log(lam)-0.5*x*x/(1+x*x);
  
  } else {

  // Potekhin et al. 1999

  double G, eta, eta0, beta, x, vc, gam;
  double ZZ=this->Ye()/this->Yi();
  
	x = this->x();
	eta=this->T8/(0.07832*sqrt(this->rho*1e-6)*this->Ye());
  	eta*=sqrt(1.0-this->Yn);
 
/*
	if (this->rho < 7.09e3*pow(1e-12*this->B,1.5)/this->Ye()) {
		double xr = 2.96e-5*this->rho*this->Ye()*1e12/this->B;			
		if (xr > x) x=xr;
	}
	*/
		//if (xr < sqrt(1.38e-8*this->T8/(9.11e-28*9e20))) xr=1e-10;
//		eta = 1.38d-8*this->T8/(9.11e-28*9e20*(sqrt(x*x+1.0)-1.0));
//	}
vc=x/sqrt(1+x*x);

  gam=this->gamma();

//  printf("T8=%lg, rho=%lg, x=%lg gam=%lg\n", this->T8, this->rho,x,gam);

  //if (gam < 1.0) return 1.0;

  beta=ZZ*PI*vc/137.0;
  eta0=0.19/pow(ZZ,1.0/6.0);
  
  // the easy part is to calculate G
  G=eta*(1.0+0.122*beta*beta)/sqrt(eta*eta+eta0*eta0);
  if (n==1) G+=0.0105*(1.0-1.0/ZZ)*(1+beta*vc*vc*vc)*eta/pow(eta*eta+0.0081,1.5);
  
  // next the s and w parameters
  double rtf, rd;
  rtf=1.0/(PI*137.0*beta);
  rd=(59.41/this->T8)*ZZ*x/(PI*3.0*137.0);

  // printf("%lg %lg %lg %lg %lg %lg ", this->T8, this->rho, beta, rd, gam, x);
  s=exp(-beta)*(rtf+rd*(1.0+0.06*gam)*exp(-sqrt(gam)));
  w=13.0*(1.0+beta/3.0)/rd;

  //printf("%lg %lg\n", s, w);
    
  // now calculate the exponential integrals    

  //  if (0) {

  if (w > 100.0) {

    L1=0.5*(log((1.0+s)/s)-1.0/(1.0+s));
    L2=(2.0*s+1.0)/(2.0*s+2.0)-s*log((1.0+s)/s);

  } else {

  if (s < 1e-2 && s < 1e-2/w) {
    L1=0.5*(expint(1,w)+log(w)+0.5772);
    L2=(exp(-w)-1.0+w)/(2.0*w);
    
  } else {
//    double I1,I2;
      // printf("s=%lg w=%lg    rho=%lg T8=%lg\n",s,w,this->rho, this->T8);
    I1=expint(1,s*w);
    //printf("I1=%lg\n",I1);
    I2=expint(1,w*(1.0+s));
    //printf("I2=%lg\n",I2);
 
    L1=log((1.0+s)/s)+(s/(1.0+s))*(1.0-exp(-w))-(1.0+s*w)*exp(s*w)*
      (I1-I2);
    L1/=2.0;
    
    L2=((exp(-w)-1.0+w)/w)-(s*s/(1.0+s))*(1.0-exp(-w))-2*s*log((1.0+s)/s)+
      s*(2.0+s*w)*exp(s*w)*(I1-I2);
    L2/=2.0;
    
    // printf("L1=%lg L2=%lg I1-I2=%lg\n", L1,L2,I1-I2);
  }
  }
  // put it all together
  G*=(L1-vc*vc*L2);

  double D=exp(-0.42*sqrt(x/((this->A[1]*(1.0-this->Yn))*this->Z[1]))*3.0*exp(-9.1*eta));
  G*=D;

	
	
	// magnetic field part
	double corr=1.0;
	if (0) {
//	if (this-> B > 0.0) {
	double bn=this->B/4.414e13;
	double xr;
	if (this->rho < 7.09e3*pow(1e-12*this->B,1.5)/this->Ye()) {
		xr = 2.96e-5*this->rho*this->Ye()*1e12/this->B;			
		//if (xr < sqrt(1.38e-8*this->T8/(9.11e-28*9e20))) xr=1e-10;
		
	} else {
	
		xr = this->x();
	}
	double gr=sqrt(1.0+xr*xr);	
	double v0 = xr/gr;
	double nu = xr*xr/(2.0*bn);
	int nmax = (int) nu;
	//double gamma = 22.75*pow(this->Z[1],2.0)*pow(1e-6*this->rho/this->A[1],1.0/3.0)/(100.0*this->T8);
	double aDW = 13.0*4.0*pow(9.0*PI*this->Z[1]/4.0,2.0/3.0)/(3.0*gam);
	double beta = this->Z[1]*PI*xr/(137.0*gr);
	double as = exp(-2.0*beta)*(13.0*(1.0+0.06*gam)*exp(-sqrt(gam))/aDW + gr/(137.0*PI*xr));
	double at = pow(sqrt(as)+pow(2.0+0.5*aDW,-1.0),2.0);
	double LL=log(1.0+1.0/at);
	double EE=(1.0-exp(-aDW))/aDW;
	double AA=(30.0-15.0*EE-(15.0-6.0*EE)*v0*v0)/
		(30.0-10.0*EE-(20.0-5.0*EE)*v0*v0);
	double BB=(1.5-0.5*EE+0.25*v0*v0/(1.0-2.0*v0*v0/3.0));
	double CC=(1.0-EE+0.75*v0*v0)/(1.0+v0*v0);
	double DD=1.0+0.06*LL*LL/(nmax*nmax);
	
	double xx=sqrt(2.0*(nu-nmax));	
	double term1 = 1.0+(sqrt(bn)/xr)*(AA/xx - BB*sqrt(xx) + CC*(xx-sqrt(xx))/nmax);
	double term2 = 0.2*EE + 0.07 + (3.0*xx*xx-1.0)/(2.0*nmax+1.5*xx*xx/pow(1.0+2.0*bn,2.0));
	
	double corr = DD/pow(term1,2.0) + pow(LL*xx*term2,2.0);
	corr = 1.0/sqrt(corr);	
//	printf("T=%lg rho =%lg corr=%lg\n",this->T8*1e8, this->rho, corr);
	}
	G *= corr;

  return G;

  }  

}  


double Eos::expint(int n, double x)
{
	if (n != 1) printf("I only know how to do Ei1(x)!\n");
	return gsl_sf_expint_E1(x);
}



double Eos::find_rho(void)
{
  double old, found, guess, rad, guess1;
  pt2Object=(void*) this;
  old=this->rho;


	if (Wrapper_find_rho_eqn(1e-6) > 0.0) return 1e-6;
 // if (2.521967e17*pow(this->T8,4)>this->P) return 1e-1;
	//printf("find_rho: %g %g\n", Wrapper_find_rho_eqn(1e-6),Wrapper_find_rho_eqn(1e15));

  found=zbrent(Wrapper_find_rho_eqn,1e-6,1e15,1e-6);
return found;
  if (0) {
  // first guess the density
  rad=2.521967e17*pow(this->T8,4);
  //  rad=0.0;
  guess=1.49e5*pow(this->P*1e-22,0.75)/this->Ye();  // NR deg electrons
  if (guess < 1e4) guess=(this->P-rad)*1.66e-24/
		     ((this->Ye()+this->Yi())*1.38e-8*this->T8);
  guess1=guess;
  while ((found = zbrent(Wrapper_find_rho_eqn,1e-1*guess,10.0*guess,1e-8))
	 <= 0.12*guess) {
    guess/=9.0;
         printf("new guess=%lg\n", guess);
	 printf("*"); fflush(stdout);
  }

  rad=2.521967e17*pow(this->T8,4);
  printf("%lg %lg %lg %lg %lg %lg\n", guess, found, rad, this->P, guess1, this->T8);
  
  }
  this->rho=old;

  if (found > 0.0) return found;
  else {
    printf("found zero density!\n");
    return 1e-1;
  }
}

double Eos::Wrapper_find_rho_eqn(double r)
{
  Eos* mySelf = (Eos*) pt2Object;
  // call member
  return mySelf->find_rho_eqn(r);
}

double Eos::find_rho_eqn(double r)
{
  this->rho=r;
  return this->ptot()-this->P;
}

double Eos::x(void)
{
  double x; 
  x=pow(this->Chabrier_EF()/511.0,2)-1.0; if (x<0.0) x=1e-10;   
  x=sqrt(x);
  return x;
}

double Eos::eta(void)
{
  return (this->Chabrier_EF()-511.0)/(8.625*this->T8);
}

double Eos::gamma(void)
{ 
  return 0.11*(this->YZ2()/this->Yi())*pow(this->rho*1e-5*Yi(),1.0/3.0)/this->T8;
  //  return 0.11*this->Z[1]*this->Z[1]*pow(this->rho*1e-5*Yi(),1.0/3.0)/this->T8;
}


double Eos::K_cond(double ef)
// Calculates the conductivity due to electron-ion and electron-electron collisions
// ef is the Fermi energy in keV
{
	// set up parameters
	double rY=this->rho*this->Ye();
	double x=this->x();
	double x2=sqrt(1+x*x); 
	double beta=x/x2;
	double gam=this->gamma();

	// Coulomb logarithm from Yakovlev and Urpin
	//double lam=log(pow(2*PI*Ye()/(3*Yi()),1.0/3.0)*sqrt(1.5+3.0/gam));
	//lam-=0.5*beta*beta;
	//this->lambda2=lam;

	// Coulomb logarithm for electron-ion collisions from Potekhin et al. 1999
	if (this->rho > 1e3)  this->lambda2=this->lamei(1); else this->lambda2=1.0;

  	// electron-electron collisions
  	// Note that Potekhin et al 1997 which is where we get the J(x,y) function has a misprint in the prefactor for f_ee
	// The correct expression is in Timmes 1992 or in Potekhin et al. (1999).
    double y=5.771e-3*sqrt(rY/x2)/this->T8;
  	f_ee=5.11e15*this->T8*this->T8*pow(x,1.5)*J(x,y)/pow(1+x*x,1.25);

	double f_c;
	
	if (gam < this->gamma_melt || this->Qimp == 900.0) { // if Q=900 treat as liquid 

    	// The electron-ion collision frequency
    	f_ei=1.76e16*this->lambda2*x2*YZ2()/Ye();
    	// The collision frequencies add
    	f_c=f_ee+f_ei;

  	} else { // solid --- NB assumes A=2Z and single species

		// The electron-ion collision frequency (ie. phonons) as given by Potekhin et al. 1999
		f_ep=1.76e16*this->lambda2*x2*YZ2()/Ye();
		// add exponential suppression when the Umklapp scatterings freeze out
		{
			double TU=2.2e8*sqrt(1e-12*this->rho)*this->Ye()*pow(this->Z[1]/60.0,1.0/3.0);
			if (this->T8<1e-8*TU) f_ep*=exp(-1e-8*TU/this->T8);
		}

		/* old phonons from Yakovlev & Urpin
		theta=0.56*sqrt(1e-9*this->rho)/this->T8;
		lam=(2-beta*beta)/(beta*sqrt(1+pow(theta/3.5,2.0)));
		lam+=pow(theta/5.1,2.0)*(3*this->lambda2-1+0.5*beta*beta)/
		(beta*pow(1+pow(theta/4.2,2.0),1.5));
		f_c=1.24e18*this->T8*lam;
		*/
		/* and from Baiko & Yakovlev 
		f_ep=9.55e16*this->T8*this->Fep(1)/beta; // phonons
		*/

		// Impurity scattering, Coulomb log from Itoh & Kohyama 1993
		double lam;
		{
			double ka, sm1;
			// the following is eq.(20) of IK93 for ka
			//ka=1.92*pow(this->Ye()/this->Yi(),1.0/3.0);
			// instead, we use the substitution  ka -> 2k/kTF and kTF is given by eq.(3) of Potekhin et al. 1999
			ka=sqrt(137.0*PI*beta);
			// and then put into eqs (10,16,17) of IK93
			sm1=0.5*log(1.0+0.4*ka*ka);
			lam=sm1*(1.0+2.5*beta*beta/(ka*ka))-0.5*beta*beta;
		}      

		 if (this->Yn > 0.0) f_eQ=1.76e16*this->Qimp*lam*x2/this->Z[1];
		 else f_eQ=1.76e16*this->Qimp*lam*x2/this->Z[1];

    	// sum of phonons and impurities and electrons
    	f_c=f_eQ+f_ep+f_ee;
  	}

	// the conductivity is then as given by Yakovlev & Urpin
	return 4.116e27*this->T8*rY/(x2*f_c);
}



double Eos::J(double x,double y)
{
	// from Potekhin, Chabrier, & Yakovlev 1997
	double x2=x*x;
	double b2=x2/(1.+x2);
	double y3=y*y*y;
	double y4=y3*y;
	return (1.+0.4*(3.+1./x2)/x2)*(y3*pow(1.+0.07414*y,-3.)*log((2.810-0.810*b2+y)/y)/3.+ pow(PI,5)*y4*pow(13.91+y,-4.)/6);
}





void Eos::potek_eos(double *P_out, double *cv_out_e, double *cv_out_i)
{
/*
	"""Magnetic EOS from Potekhin & Chabrier 2013. Input: Z,A,rho,T,B,LIQSOL
	where
	if LIQSOL=0 or 1 on the input, then:
	*        if (either GAMI or 1/RS is below its critical value) then
	*           liquid regime and LIQSOL=0 on the output
	*        otherwise solid regime and LIQSOL=1 on the output
	*     if LIQSOL=2 or 3 on the input, then:
	*        if (LIQSOL=2) then liquid regime and LIQSOL=2 on the output
	*        if (LIQSOL=3) then solid regime and LIQSOL=3 on the output
	*     if LIQSOL=4 or 5 on the input, then
	*         consider non-ideal free energy FC1:
	*        if (FC1=min in liquid) then liquid regime, output LIQSOL=4
	*        if (FC1=min in solid) then solid regime, output LIQSOL=5
	Output:
	*         DENS - electron number density [in a.u.=6.7483346e24 cm^{-3}]
	*         GAMI - ion-ion Coulomb coupling constant
	*         CHI = mu_e/kT, where mu_e is the electron chem.potential
	*         TPT - ionic quantum parameter (T_p/T)
	*         LIQSOL - regulator of the solid vs liquid regime (see below)
	*         PnkT - pressure / n_i kT, where n_i is the ion number density
	*         UNkT - internal energy per kT per ion
	*         SNk - total dimensionless entropy per 1 ion (assumes spin=1/2)
	*         CV - heat capacity per ion, divided by the Boltzmann constant
	*         CHIR - inverse compressibility -(d ln P / d ln V)_T ("\chi_r")
	*         CHIT = (d ln P / d ln T)_V ("\chi_T")	
	"""
*/	
	double Zion = this->Z[1];
	double CMI = this->A[1];
	double RR = this->rho;
	double TT = this->T8*1e8/3.1577e5;
	//double BB = this->B/2.3505e9;
	double GAMAG = 0.0;
	double DENS, GAMI, CCHI, TPT, LIQSOL=1, PnkT, UNkT,SNk,CCVI,CCVE,CHIR,CHIT;
//	if (TT<0.0) TT=100.0;
	eosmag_(&Zion,&CMI,&RR,&TT,&GAMAG,&DENS,&GAMI,&CCHI,&TPT,&LIQSOL,&PnkT,&UNkT,&SNk,&CCVE,&CCVI,
			&CHIR,&CHIT);
	//Multiply pressure by 8.31447e13 rho T6/CMImean to get cgs pressure
	*P_out = PnkT * 8.31447e13 * this->rho * 100.0*this->T8/this->A[1];
	*cv_out_e = CCVE * 8.31447e7/this->A[1];
	*cv_out_i = CCVI * 8.31447e7/this->A[1];
}


double Eos::potek_cond(void)
// returns the thermal conductivity in cgs from Potekhin's fortran code
{
	double s1,s2,s3,k1,k2,k3;
	//double null=0.0;
	double Zimp=sqrt(this->Qimp), AA=this->A[1]*(1.0-this->Yn);
	double Bfield=this->B/4.414e13;
	double temp=this->T8*1e2/5930.0;
	double rr=this->rho/(this->A[1]*15819.4*1822.9);
	condegin_(&temp,&rr,&Bfield,&this->Z[1],&AA,&this->A[1],&Zimp, &s1,&s2,&s3,&k1,&k2,&k3);
	this->Kperp = k2*2.778e15;
	return k1*2.778e15;

// This was my attempt at SF phonons:
/*		double ksph =0.0;
		if (EOS.Yn > 0.0 && G.include_sph) {
			// now SF conductivity
			double lsph = 0.001 * pow(EOS.rho/4e11,-1.0/3.0);  // mfp in cm

			double vs = 1.05e-27/(sqrt(3.0)*1.67e-24*3e10);
			vs*=pow(3.0*PI*PI*EOS.rho*EOS.Yn/1.67e-24,1.0/3.0);	
			//printf("%lg %lg %lg %lg\n", EOS.T8, EOS.rho, EOS.Yn, vs);

	//		double vs = 0.1 * pow(EOS.rho/4e11,1.0/3.0);
			ksph = 1.5e22 * pow(EOS.T8,3.0) * pow(0.1/vs,2.0) * lsph;
		}

		//ksph=0.0;

		*Kperp += ksph;
		return kk + ksph;
		*/
}





double Eos::econd(void)
  // calculates the electrical conductivity
{
  double x1, x2, sig, x, lambda, nu, beta;

  if (this->gamma() < this->gamma_melt || this->Qimp == 900.0) { // if Q=900 treat as liquid
    
    // This is the method from the WD paper, where I interpolate using x
    // choose appropriate value for x
    x1=this->x();
    x2=0.26*sqrt(this->T8);
    x=sqrt(x1*x1+x2*x2);
    //if (x1>x2) x=x1; else x=x2;

  //  x=x1;
    sig=8.48e21*this->Ye()*pow(x,3.0)/(this->YZ2()*(1+x*x));
    //printf("got here\n");

    sig/=this->lamei(0);
    
    //    printf("back\n");

    // here, write sig directly in terms of Fermi integrals
    //    sig=9.47e23*pow(this->T8,3.0)*this->Fermi(2.0,this->eta())/this->rho;
    //sig/=this->YZ2()*this->lamei();
    
  } else { // solid --- NB assumes A=2Z and single species
    double TU,ka,sm1;

    TU=2.2e8*sqrt(1e-12*this->rho)*this->Ye()*pow(this->Z[1]/60.0,1.0/3.0);

    x=this->x(); beta=x/sqrt(1+x*x);
    nu=9.55e16*this->T8*this->Fep(0)/beta; // phonons

    // add exponential suppression when the Umklapp scatterings freeze out
    //if (this->T8  < 1e-8*TU) 
    nu*=exp(-1e-8*TU/this->T8);

    //nu=9.55e16*this->T8*13.0/beta;
    /* old phonons from Urpin & Yakovlev
    theta=0.56*sqrt(1e-9*this->rho)/this->T8;
    nu=1.24e18*this->T8*(2-beta*beta)/(beta*sqrt(1+pow(theta/3.5,2.0)));
    */

    // Coulomb log from Itoh & Kohyama 1996
    ka=1.92*pow(this->Ye()/this->Yi(),1.0/3.0);
    sm1=0.5*log(1.0+0.4*ka*ka);
    lambda=sm1*(1.0+2.5*beta*beta/(ka*ka))-0.5*beta*beta;

    nu+=1.75e16*this->Qimp*lambda*sqrt(1+x*x)/this->Z[1]; // impurities

    //sig=1.49e22*x*x*beta*1e16/nu;
    sig=1.52e25*1e17*pow(this->rho*1e-12*Ye(),2.0/3.0)/nu;
  }  
  return sig;
}




