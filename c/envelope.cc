#include "../h/envelope.h"
#include "math.h"
#include <stdlib.h>

void* pt2EnvObject;

Envelope::Envelope()
{	
	pt2EnvObject=(void*) this;
  	
	static Eos myeos(1);
    this->EOS = &myeos;
	this->EOS->X[1]=1.0;
	this->EOS->Qimp=1.0;
	this->EOS->accr=1;	

	this->use_potek_eos_in_He=1;
	this->use_potek_eos_in_Fe=0;
	this->use_potek_cond_in_He=1;
	this->use_potek_cond_in_Fe=1;
	this->use_potek_kff=0;

	this->ODE.init(1,dynamic_cast<Ode_Int_Delegate *>(this));
	this->ODE2.init(1,dynamic_cast<Ode_Int_Delegate *>(this));
		
	// Default values
	this->g = 2.28e14;
	
}


Envelope::~Envelope() 
{
	this->ODE.tidy();
	this->ODE2.tidy();
}


void Envelope::make_grid(double yi, double B)
// yi is the log10 of the base column of the He layer
// B is the magnetic field strength
{
	this->fp = fopen("out/grid","w");
	this->yi = yi;
	this->Bfield = B;
	this->EOS->B=this->Bfield;

	for (int i=0; i<=90; i++) {
		double flux=17.0+i*0.1;
		this->F=pow(10.0,flux);
		this->doint();
		printf("."); fflush(stdout);
     	//for (int j=1; j<=ODE2.kount; j++) {
		//	fprintf(this->fp, "%lg %lg %lg %lg\n", this->ODE2.get_x(j), log10(this->ODE2.get_y(1,j)), flux, this->yt);
		//}
	   	for (int j=1; j<=ODE.kount; j++) {
				fprintf(this->fp, "%lg %lg %lg %lg\n", this->ODE.get_x(j), log10(this->ODE.get_y(1,j)), flux, this->yt);
		}
	}

	printf("\n");
	fclose(this->fp);
}


void Envelope::doint(void)
// for the specified flux, integrate inwards to see if we match the base temperature
{
  	// we do this in two steps: light element layer first
	this->EOS->use_potek_eos = use_potek_eos_in_He;
	this->EOS->use_potek_cond = use_potek_cond_in_He;
	this->EOS->use_potek_kff = use_potek_kff;
	double yi;
	if (this->yi == 0.0) {
		this->EOS->A[1]=56.0; this->EOS->Z[1]=26.0; // Iron		
		yi = 8.0;
	} else {
		this->EOS->A[1]=4.0; this->EOS->Z[1]=2.0; // Helium
		yi = this->yi;
	}
  	// set surface temperature. We integrate from tau=2/3
	double y1=1e-3,y2=1e8;
	this->yt=zbrent(this->Wrapper_find_surf_eqn,y1,y2,1e-6);
  	if (this->yt==y1 || this->yt==y2) printf("yt out of bounds (%lg)\n", this->yt);
	if (this->yt > pow(10.0,yi)) {
		printf("The photosphere lies beneath the helium column! Increase yi and try again.\n");
		exit(0);
	}

  	double Tt=pow(this->F/5.67e-5,0.25);
	this->ODE2.set_bc(1,Tt);

  	// integrate
  	this->ODE2.go_simple(log10(this->yt),yi,(int)((yi-log10(this->yt))/0.1));

  	// keep the base temperature for the next integration
  	double base_T=this->ODE2.get_y(1,this->ODE2.kount);

  	// setup ocean
  	this->ODE.set_bc(1,base_T);
	this->EOS->use_potek_cond = use_potek_cond_in_Fe;
	this->EOS->use_potek_eos = use_potek_eos_in_Fe;
 	this->EOS->A[1]=56.0; this->EOS->Z[1]=26.0;  // Iron
  
  	// integrate through the ocean to the desired depth 
	this->ODE.go_simple(yi,18.5,(int)((18.5-yi)/0.1));
	
	// reset EOS for next time
	this->EOS->Yn=0.0; this->EOS->set_Ye=0.0; 
}



double Envelope::Wrapper_find_surf_eqn(double y)
{
  Envelope* mySelf = (Envelope*) pt2EnvObject;
  return mySelf->find_surf_eqn(y);
}

double Envelope::find_surf_eqn(double y)
{
	this->EOS->P=this->g*y; this->EOS->T8=1e-8*pow(this->F/5.67e-5,0.25);
	this->EOS->rho=this->EOS->find_rho();
	//printf("%lg %lg %lg\n",this->EOS->P,this->EOS->rho,this->EOS->T8);
	return this->EOS->opac()*y-2.0/3.0;
}


void Envelope::derivs(double x, double ff[], double dfdx[])
// Evaluate derivative dT/dy for the integrator
{
  	double y=pow(10.0,x);
  	double T=ff[1];

  	// find density; the composition is already set
  	this->EOS->P=this->g*y; EOS->T8=T*1e-8; 
	if (x > this->yi) EOS->set_composition_by_pressure();
  	this->EOS->rho=EOS->find_rho();

	double kappa=EOS->opac();

  	// Heat equation for constant flux
  	dfdx[1]=2.303*y*3305.1*this->F*kappa/pow(T,3.0);

	// if the temperature gradient is superadiabatic, set it to the adiabatic gradient (convection)
	if (this->EOS->gamma() > this->EOS->gamma_melt) {
		double conv_grad = 2.303*T*EOS->del_ad();
		if (dfdx[1] > conv_grad) dfdx[1]=conv_grad;
	}
}

void Envelope::jacobn(double, double *, double *, double **, int) 
{
}

