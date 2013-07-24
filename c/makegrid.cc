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
double find_surf_eqn(double y);

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
	//for (int i=2; i<=5; i++) {
	//	{ int i = 2;
	  //  	double F=17.0+i*1.2;
	for (int i=0; i<=180; i++) {
	    	double F=17.0+i*0.05;
   		doint(F);
		printf("Did doint for F=%lg\n",F);
    	printf("."); fflush(stdout);

		EOS.init(1);
		EOS.B=G.Bfield;
		EOS.use_potek_cond = 0;
		EOS.use_potek_eos = 0;

	  	EOS.X[1]=1.0; EOS.A[1]=4.0; EOS.Z[1]=2.0;

    	for (int j=1; j<=ODE2.kount; j++) {


					double y=pow(10.0,ODE2.get_x(j));
				  	double T=ODE2.get_y(1,j);

			  	// find density; the composition is already set
			  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
			  	EOS.rho=EOS.find_rho();
				double rhorho=EOS.rho;

		      	printf("%lg %lg %lg %lg %lg\n", ODE2.get_x(j), log10(ODE2.get_y(1,j)), F, log10(rhorho),
						log10(EOS.eps_nu()));

	      	fprintf(fp.out, "%lg %lg %lg %lg %lg\n", ODE2.get_x(j), log10(ODE2.get_y(1,j)), F, log10(rhorho),
					log10(EOS.eps_nu()));
	} 	
	
	EOS.use_potek_cond = 1;
	EOS.use_potek_eos = 1;
	EOS.A[1]=56.0; EOS.Z[1]=26.0;  EOS.X[1]=1.0;
	

    	for (int j=1; j<=ODE.kount; j++) {
			double y=pow(10.0,ODE.get_x(j));
		  	double T=ODE.get_y(1,j);

		  	// find density; the composition is already set
		  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
		  	EOS.rho=EOS.find_rho();
			if (log10(y) > G.yi) EOS.set_comp();
			double rhorho=EOS.rho;

	    	printf("%lg %lg %lg %lg %lg\n", ODE.get_x(j), log10(ODE.get_y(1,j)), F, log10(rhorho),
					log10(EOS.eps_nu()*y));

      	fprintf(fp.out, "%lg %lg %lg %lg %lg\n", ODE.get_x(j), log10(ODE.get_y(1,j)), F, log10(rhorho),
				log10(EOS.eps_nu()*y));
}
	EOS.tidy();
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
	EOS.use_potek_eos = 0;
	EOS.use_potek_cond = 0;
	
  	EOS.X[1]=1.0; EOS.A[1]=4.0; EOS.Z[1]=2.0;
	//EOS.X[1]=1.0; EOS.A[1]=56.0; EOS.Z[1]=26.0;

  	// set surface temperature. We integrate from tau=2/3
  	double Tt,yt;
  	yt=zbrent(find_surf_eqn,1e-5,1e5,1e-5);
  	if (yt==1e-5 || yt==1e3) printf("yt out of bounds (%lg)\n", yt);
	//printf("yt=%lg\n", yt);
  	Tt=pow(F/5.67e-5,0.25);
  	ODE2.set_bc(1,Tt);

	printf("Starting heliun layer\n");


  	// integrate
//  	ODE2.go(log10(yt),G.yi,1e-6,1e-8,lc_derivs);
  	ODE2.go_simple(log10(yt),G.yi,(int)(80*(18.5-G.yi)),lc_derivs);

  	// keep the base temperature for the next integration
  	base_T=ODE2.get_y(1,ODE2.kount);

	printf("Finished heliun layer\n");

  	// tidy up and reinitialize for the ocean
  	ODE.set_bc(1,base_T);
  	EOS.tidy();
  	EOS.init(1);
	EOS.B=G.Bfield;
	EOS.Q=0.0;
	EOS.use_potek_cond = 1;
	EOS.use_potek_eos = 1;
 	EOS.A[1]=56.0; EOS.Z[1]=26.0;  EOS.X[1]=1.0;
  
  	// integrate through the ocean to the desired depth 
  	ODE.go_simple(G.yi,18.5,(int)(80*(18.5-G.yi)),lc_derivs);

/*
	int flag =0;
	for (int ii=1; ii<=ODE.kount; ii++) {
		double y=pow(10.0,ODE.get_x(ii));
	  	double T=ODE.get_y(1,ii);

	  	// find density; the composition is already set
	  	EOS.P=G.g*y; EOS.T8=T*1e-8; 
	  	EOS.rho=EOS.find_rho();
		if (log10(y) > G.yi) EOS.set_comp();

	  //	double kappa=EOS.opac();
	
		//if (EOS.kcond < EOS.kappa_rad && !flag) {
		//	flag=1;
			//printf("%lg %lg %lg %lg %lg %lg\n", y, EOS.rho, EOS.T8*1e8, EOS.kappa_rad, EOS.kcond,3.024e20*pow(EOS.T8,3)/(potek_cond()*EOS.rho));
//		}

	  	// here I call Potekhin's routine for the conductivity  
	 	//EOS.kcond=3.024e20*pow(EOS.T8,3)/(potek_cond()*EOS.rho); 
	  	//kappa=1.0/((1.0/EOS.kcond)+(1.0/EOS.kappa_rad));
	 	
	}
*/

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
	EOS.Yn=0.0;
	if (x > G.yi) EOS.set_comp();

  	double kappa=EOS.opac();
	if (0) {//}(EOS.kcond < EOS.kappa_rad*5.0) {
		double mykcond = EOS.kcond;
		// here I call Potekhin's routine for the conductivity  
 		EOS.kcond=3.024e20*pow(EOS.T8,3)/(EOS.potek_cond()*EOS.rho);
		// and then swap it in
  		kappa=1.0/((1.0/EOS.kcond)+(1.0/EOS.kappa_rad));
		//printf("%lg %lg %lg %lg %lg %lg\n", EOS.rho, EOS.T8, EOS.kcond, mykcond, EOS.kcond/mykcond, potek_cond());
	}

  	// Heat equation for constant flux
  	dfdx[1]=2.303*y*3305.1*G.F*kappa/pow(T,3.0);

//	if (EOS.B>0.0) {
//		double conv_grad = 2.303*T*EOS.del_ad();
//		if (dfdx[1] > conv_grad) dfdx[1]=conv_grad;
//	}

// 	printf("%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", x, ff[1], dfdx[1], 
//	EOS.opac(), EOS.Chabrier_EF(), EOS.kff, EOS.kes, EOS.kcond, log10(EOS.rho), 
//		EOS.kgaunt, EOS.x(),EOS.del_ad());
}

