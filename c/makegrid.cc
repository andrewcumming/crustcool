// makegrid.cc
//
// The idea here is to calculate a range of constant flux solutions 
// on a constant grid in y, so that we get a mapping between 
// T_b,y_b and the flux
 
#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>

#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/odeint.h"
#include "../h/eos.h"
#include "../h/spline.h"
#include "../h/useful.h"

// -------- Global variables ----------

Eos EOS;
Ode_Int ODE;
Ode_Int ODE2;

struct global {
  double g;
  double F;
double Bfield;
  double yi; 
} G;

struct file_pointers {
  FILE *out;
} fp;


// --------- Declarations ------------
void jacobn(double, double *, double *, double **, int){};
void doint(double F);
void lc_derivs(double y, double ff[], double dfdy[]);
extern "C"{
  void condegin_(double *temp,double *densi,double *B,double *Zion,double *CMI,double *CMI1,double *Zimp, double *RSIGMA,double *RTSIGMA,double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA);
 // void conduct_(double *temp,double *densi,double *B,double *Zion,double *CMI,double *Zimp, double *RSIGMA,double *RTSIGMA,double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA,double *RHKAPPA);
}
double find_surf_eqn(double y);
double potek_cond(void);

//------------------------------------------------------------------------


int main(void)
{
  	fp.out=fopen("out/grid","w");

  	G.g=2.28e14;  // choose gravity to match Ed, we can rescale to arbitrary g later

  	printf("Enter log10 base column of He layer ..."); scanf("%lg",&G.yi);

  	EOS.init(1);
	EOS.Q=0.0;  EOS.accr=1; // set crust composition to accreted crust
  	EOS.X[1]=1.0; 
 
	printf("Enter B field in G (0 for unmagnetized)..."); scanf("%lg",&G.Bfield);


  	ODE.init(1);
  	ODE2.init(1);

  	// now try different surface fluxes
 //	for (int i=0; i<=5; i++) {
   // 	double F=17.0+i*1.25;
	for (int i=0; i<=125; i++) {
	    	double F=17.0+i*0.05;
   		doint(F);
    	printf("."); fflush(stdout);
    	//  for (int j=1; j<=ODE2.kount; j++) 
    	// fprintf(fp.out, "%lg %lg %lg\n", ODE2.get_x(j), log10(ODE2.get_y(1,j)), F);
	 	
	
	//EOS.init(1);
	//	EOS.B=G.Bfield;
	  //	EOS.A[1]=56.0; EOS.Z[1]=26.0;  EOS.X[1]=1.0;

    	for (int j=1; j<=ODE.kount; j++) {
			double y=pow(10.0,ODE.get_x(j));
		  	double T=ODE.get_y(1,j);

		  	// find density; the composition is already set
		//  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
		  //	EOS.rho=EOS.find_rho();
			//if (log10(y) > G.yi) EOS.set_comp();
//			double rhorho=EOS.rho;
			double rhorho=0.0;

      	fprintf(fp.out, "%lg %lg %lg %lg\n", ODE.get_x(j), log10(ODE.get_y(1,j)), F, rhorho);
}
//EOS.tidy();
  	}
  	printf("\n");

  	// tidy up
  	fclose(fp.out);
  	ODE.tidy();
  	ODE2.tidy();
}


double find_surf_eqn(double y)
{
	EOS.P=G.g*y; EOS.T8=1e-8*pow(G.F/5.67e-5,0.25);
	EOS.rho=EOS.find_rho();
	return EOS.opac()*y-2.0/3.0;
}


void doint(double F) 
{
  	// for a particular choice of flux, integrate inwards to see
  	// if we match the base temperature

  	F=pow(10.0,F);   // user gives log_10 F
  	G.F=F;
  	double base_T;

  	// we do this in two steps: light element layer
  	EOS.init(1); 
	EOS.B=G.Bfield;
  	EOS.X[1]=1.0; EOS.A[1]=4.0; EOS.Z[1]=2.0;
	//EOS.X[1]=1.0; EOS.A[1]=56.0; EOS.Z[1]=26.0;

  	// set surface temperature. We integrate from tau=2/3
  	double Tt,yt;
  	yt=zbrent(find_surf_eqn,1e-3,1e2,1e-5);
  	if (yt==1e-3 || yt==1e2) printf("yt out of bounds (%lg)\n", yt);
	//printf("yt=%lg\n", yt);
  	Tt=pow(F/5.67e-5,0.25);
  	ODE2.set_bc(1,Tt);

  	// integrate
  	ODE2.go(log10(yt),G.yi,1e-6,1e-8,lc_derivs);

  	// keep the base temperature for the next integration
  	base_T=ODE2.get_y(1,ODE2.kount);

  	// tidy up and reinitialize for the ocean
  	ODE.set_bc(1,base_T);
  	EOS.tidy();
  	EOS.init(1);
	EOS.B=G.Bfield;
  	EOS.A[1]=56.0; EOS.Z[1]=26.0;  EOS.X[1]=1.0;
  
  	// integrate through the ocean to the desired depth 
  	ODE.go_simple(G.yi,18.5,(int)(40*(18.5-G.yi)),lc_derivs);

	int flag =0;
	for (int ii=1; ii<=ODE.kount; ii++) {
		double y=pow(10.0,ODE.get_x(ii));
	  	double T=ODE.get_y(1,ii);

	  	// find density; the composition is already set
	  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
	  	EOS.rho=EOS.find_rho();
		if (log10(y) > G.yi) EOS.set_comp();

	  	double kappa=EOS.opac();
	
		if (EOS.kcond < EOS.kappa_rad && !flag) {
			flag=1;
			//printf("%lg %lg %lg %lg %lg %lg\n", y, EOS.rho, EOS.T8*1e8, EOS.kappa_rad, EOS.kcond,3.024e20*pow(EOS.T8,3)/(potek_cond()*EOS.rho));
		}

	  	// here I call Potekhin's routine for the conductivity  
	 	//EOS.kcond=3.024e20*pow(EOS.T8,3)/(potek_cond()*EOS.rho); 
	  	//kappa=1.0/((1.0/EOS.kcond)+(1.0/EOS.kappa_rad));
	 	
	}


  	EOS.tidy();
}


void lc_derivs(double x, double ff[], double dfdx[])
// Evaluate derivative dT/dy
{
  	double y=pow(10.0,x);
  	double T=ff[1];

  	// find density; the composition is already set
  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
  	EOS.rho=EOS.find_rho();
	if (x > G.yi) EOS.set_comp();

  	double kappa=EOS.opac();

  	// here I call Potekhin's routine for the conductivity  
 	//EOS.kcond=3.024e20*pow(EOS.T8,3)/(potek_cond()*EOS.rho); 
  	//kappa=1.0/((1.0/EOS.kcond)+(1.0/EOS.kappa_rad));

  	// Heat equation for constant flux
  	dfdx[1]=2.303*y*3305.1*G.F*kappa/pow(T,3.0);

  	//printf("%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", x, ff[1], dfdx[1], 
	//EOS.opac(), EOS.Chabrier_EF(), EOS.kff, EOS.kes, EOS.kcond, log10(EOS.rho), 
	//EOS.kgaunt, EOS.x());
}


double potek_cond()
// returns the thermal conductivity in cgs from Potekhin's fortran code
{
      double s1,s2,s3,k1,k2,k3;
      double null=0.0, Zimp=sqrt(EOS.Q), AA=EOS.A[1]*(1.0-EOS.Yn);
	  double Bfield=EOS.B/4.414e13;
      double temp=EOS.T8*1e2/5930.0;
      double rr=EOS.rho/(EOS.A[1]*15819.4*1822.9);
   //   conduct_(&temp,&rr,&Bfield,&EOS.Z[1],&AA,&Zimp, &s1,&s2,&s3,&null,&k1,&k2,&k3,
	//			&null,&null,&null,&null,&null,&null,&null,&null);
  	condegin_(&temp,&rr,&Bfield,&EOS.Z[1],&AA,&EOS.A[1],&Zimp, &s1,&s2,&s3,&k1,&k2,&k3);
      return k1*2.778e15;
}
